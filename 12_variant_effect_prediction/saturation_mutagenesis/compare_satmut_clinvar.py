"""
Filters the 31.5M-mutation saturation scan down to the subset that
matches a real observed ClinVar missense variant (matched by
gene_symbol + position + wt_aa + mut_aa -- ClinVar's RefSeq-transcript
numbering is assumed to line up with our UniProt canonical-sequence
numbering; a real caveat, not verified per-gene), then compares
TANGO/protGPS/DeepLoc change-from-WT between:
  - ClinVar Pathogenic/Likely pathogenic missense matches
  - ClinVar Benign/Likely benign missense matches
  - ClinVar Uncertain-significance / other missense matches
  - A background sample of the full unfiltered 31.5M-mutation set
    (every mutation, not just ones ClinVar happens to have observed)

Outputs per-row delta values (not just means) for each bucket so the
comparison plots can show real distributions, not just averages.
"""
import csv
import json
import re
import sys

csv.field_size_limit(sys.maxsize)

BASE = "/u/project/kappel/ddcohn/protein_variant_effects/rbp_only"
CLINVAR_CHANGES = "/u/project/kappel/ddcohn/protein_variant_effects/clinvar_protein_changes.tsv"
CLINVAR_META_RBP = f"{BASE}/clinvar_metadata_lookup_rbp.tsv"
SATMUT_ALL = f"{BASE}/satmut_all_predictions.tsv"
WT_TANGO = f"{BASE}/wt_tango_scores.tsv"
MERGED_TABLE = "/u/project/kappel/ddcohn/RBP/merged_260903B_ptm_pred_slim_RBPs_domains_spurs_rsa_iso_hpa_with_rnp_rna_class.csv"
OUT_JSON = f"{BASE}/satmut_vs_clinvar.json"

PROTGPS_COLS = [
    "NUCLEAR_SPECKLE_Score", "P-BODY_Score", "PML-BDOY_Score",
    "POST_SYNAPTIC_DENSITY_Score", "STRESS_GRANULE_Score", "CHROMOSOME_Score",
    "NUCLEOLUS_Score", "NUCLEAR_PORE_COMPLEX_Score", "CAJAL_BODY_Score",
    "RNA_GRANULE_Score", "CELL_JUNCTION_Score", "TRANSCRIPTIONAL_Score",
]
DEEPLOC_NUMERIC_COLS = [
    "Cytoplasm", "Nucleus", "Extracellular", "Cell membrane", "Mitochondrion",
    "Plastid", "Endoplasmic reticulum", "Lysosome/Vacuole", "Golgi apparatus",
    "Peroxisome", "Peripheral", "Transmembrane", "Lipid anchor", "Soluble",
]

AA3_TO_1 = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q",
    "Glu": "E", "Gly": "G", "His": "H", "Ile": "I", "Leu": "L", "Lys": "K",
    "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W",
    "Tyr": "Y", "Val": "V",
}

MUT_ID_RE = re.compile(r"^(?P<acc>[^|]+)\|(?P<gene>[^|]+)\|(?P<wt>[A-Z])(?P<pos>\d+)(?P<mut>[A-Z])$")

PATHOGENIC = {"Pathogenic", "Likely pathogenic", "Pathogenic/Likely pathogenic"}
BENIGN = {"Benign", "Likely benign", "Benign/Likely benign"}


def bucket_significance(sig):
    if sig in PATHOGENIC:
        return "pathogenic"
    if sig in BENIGN:
        return "benign"
    return "uncertain"


def load_clinvar_keys():
    var_sig = {}
    with open(CLINVAR_META_RBP, newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            var_sig[row["VariationID"]] = row["ClinicalSignificance"]
    print(f"Loaded {len(var_sig)} RBP ClinVar VariationID->significance entries")

    keys = {}  # (gene, pos, wt, mut) -> bucket
    n_missense = 0
    n_matched_rbp = 0
    with open(CLINVAR_CHANGES, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            if row["category"] != "missense":
                continue
            n_missense += 1
            vid = row["VariationID"]
            sig = var_sig.get(vid)
            if sig is None:
                continue  # not an RBP-gene variant
            n_matched_rbp += 1
            try:
                details = json.loads(row["details"])
                pos = details["pos"]
                wt = AA3_TO_1.get(details["wt"])
                mut = AA3_TO_1.get(details["mut"])
            except (json.JSONDecodeError, KeyError):
                continue
            if wt is None or mut is None:
                continue
            gene = row["GeneSymbol"]
            key = (gene, pos, wt, mut)
            bucket = bucket_significance(sig)
            keys[key] = bucket

    print(f"Total missense ClinVar variants: {n_missense}")
    print(f"Missense variants matched to an RBP gene: {n_matched_rbp}")
    print(f"Unique (gene,pos,wt,mut) ClinVar missense keys usable for matching: {len(keys)}")
    from collections import Counter
    print(f"Bucket breakdown: {Counter(keys.values())}")
    return keys


def load_wt_tango():
    wt = {}
    with open(WT_TANGO, newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            try:
                wt[row["uniprot_accession"]] = float(row["AGG"])
            except ValueError:
                pass
    return wt


def load_wt_protgps_deeploc():
    protgps_wt = {}
    deeploc_wt = {}
    seen = set()
    with open(MERGED_TABLE, newline="") as f:
        for row in csv.DictReader(f):
            acc = row["uniprot_accession"]
            if acc in seen:
                continue
            seen.add(acc)
            pvals = []
            ok = True
            for c in PROTGPS_COLS:
                v = row.get(f"protgps_{c}", "").strip()
                if not v or v == "NA":
                    ok = False
                    break
                pvals.append(float(v))
            if ok:
                protgps_wt[acc] = pvals
            dvals = []
            ok = True
            for c in DEEPLOC_NUMERIC_COLS:
                v = row.get(f"deeploc_{c}", "").strip()
                if not v or v == "NA":
                    ok = False
                    break
                dvals.append(float(v))
            if ok:
                deeploc_wt[acc] = dvals
    return protgps_wt, deeploc_wt


def main():
    clinvar_keys = load_clinvar_keys()
    wt_tango = load_wt_tango()
    wt_protgps, wt_deeploc = load_wt_protgps_deeploc()
    print(f"WT baselines: {len(wt_tango)} TANGO, {len(wt_protgps)} protGPS, {len(wt_deeploc)} DeepLoc")

    buckets = {"pathogenic": {"tango": [], "protgps": [], "deeploc": []},
               "benign": {"tango": [], "protgps": [], "deeploc": []},
               "uncertain": {"tango": [], "protgps": [], "deeploc": []},
               "background": {"tango": [], "protgps": [], "deeploc": []}}

    SAMPLE_STRIDE = 60
    n_total = 0
    n_clinvar_matched = 0

    with open(SATMUT_ALL, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for i, row in enumerate(reader):
            n_total += 1
            m = MUT_ID_RE.match(row["mutation_id"])
            if not m:
                continue
            acc = m.group("acc")
            gene = m.group("gene")
            pos = int(m.group("pos"))
            wt_aa = m.group("wt")
            mut_aa = m.group("mut")

            key = (gene, pos, wt_aa, mut_aa)
            bucket_name = clinvar_keys.get(key)

            is_background_sample = (i % SAMPLE_STRIDE == 0)
            if bucket_name is None and not is_background_sample:
                continue

            tango_d = None
            wt_agg = wt_tango.get(acc)
            if wt_agg is not None:
                try:
                    tango_d = abs(float(row["AGG"]) - wt_agg)
                except ValueError:
                    pass

            protgps_d = None
            wt_p = wt_protgps.get(acc)
            if wt_p is not None:
                try:
                    protgps_d = sum(abs(float(row[c]) - w) for c, w in zip(PROTGPS_COLS, wt_p))
                except ValueError:
                    pass

            deeploc_d = None
            wt_d = wt_deeploc.get(acc)
            if wt_d is not None:
                try:
                    deeploc_d = sum(abs(float(row[c]) - w) for c, w in zip(DEEPLOC_NUMERIC_COLS, wt_d))
                except ValueError:
                    pass

            if bucket_name is not None:
                n_clinvar_matched += 1
                if tango_d is not None:
                    buckets[bucket_name]["tango"].append(tango_d)
                if protgps_d is not None:
                    buckets[bucket_name]["protgps"].append(protgps_d)
                if deeploc_d is not None:
                    buckets[bucket_name]["deeploc"].append(deeploc_d)

            if is_background_sample:
                if tango_d is not None:
                    buckets["background"]["tango"].append(tango_d)
                if protgps_d is not None:
                    buckets["background"]["protgps"].append(protgps_d)
                if deeploc_d is not None:
                    buckets["background"]["deeploc"].append(deeploc_d)

    print(f"\nTotal satmut rows: {n_total}")
    print(f"Matched a real ClinVar missense variant: {n_clinvar_matched}")
    for b in buckets:
        print(f"  {b}: tango n={len(buckets[b]['tango'])}, protgps n={len(buckets[b]['protgps'])}, deeploc n={len(buckets[b]['deeploc'])}")

    with open(OUT_JSON, "w") as f:
        json.dump({"n_total": n_total, "n_clinvar_matched": n_clinvar_matched, "buckets": buckets}, f)
    print(f"\nWrote {OUT_JSON}")


if __name__ == "__main__":
    main()
