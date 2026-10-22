"""
Extends compare_satmut_clinvar.py to also include COSMIC/CMC: COSMIC's
raw GenomeScreensMutant download has no pathogenicity/significance tier
(that's a separate, licensed COSMIC product not present in this data),
so per user instruction, every COSMIC missense variant that matches an
RBP gene is used as one undifferentiated "observed in cancer" bucket --
no pathogenic/benign split for COSMIC, unlike ClinVar.

Buckets compared:
  - clinvar_pathogenic / clinvar_benign (real curated ClinVar labels)
  - cosmic (all COSMIC missense matches, no significance label available)
  - background (systematic sample of the full 31.5M-mutation scan)
"""
import csv
import json
import re
import sys

csv.field_size_limit(sys.maxsize)

BASE = "/u/project/kappel/ddcohn/protein_variant_effects/rbp_only"
CLINVAR_CHANGES = "/u/project/kappel/ddcohn/protein_variant_effects/clinvar_protein_changes.tsv"
CLINVAR_META_RBP = f"{BASE}/clinvar_metadata_lookup_rbp.tsv"
CMC_CHANGES = "/u/project/kappel/ddcohn/protein_variant_effects/cmc_protein_changes.tsv"
SATMUT_MANIFEST = f"{BASE}/satmut_manifest.tsv"
SATMUT_ALL = f"{BASE}/satmut_all_predictions.tsv"
WT_TANGO = f"{BASE}/wt_tango_scores.tsv"
MERGED_TABLE = "/u/project/kappel/ddcohn/RBP/merged_260903B_ptm_pred_slim_RBPs_domains_spurs_rsa_iso_hpa_with_rnp_rna_class.csv"
OUT_JSON = f"{BASE}/satmut_vs_clinvar_cosmic.json"

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


def load_rbp_genes():
    genes = set()
    with open(SATMUT_MANIFEST, newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            genes.add(row["RBP_name"])
    print(f"RBP gene set: {len(genes)}")
    return genes


def load_clinvar_keys():
    var_sig = {}
    with open(CLINVAR_META_RBP, newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            var_sig[row["VariationID"]] = row["ClinicalSignificance"]

    path_keys, benign_keys = set(), set()
    with open(CLINVAR_CHANGES, newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if row["category"] != "missense":
                continue
            sig = var_sig.get(row["VariationID"])
            if sig not in PATHOGENIC and sig not in BENIGN:
                continue
            try:
                d = json.loads(row["details"])
                wt = AA3_TO_1.get(d["wt"])
                mut = AA3_TO_1.get(d["mut"])
                pos = d["pos"]
            except (json.JSONDecodeError, KeyError):
                continue
            if wt is None or mut is None:
                continue
            key = (row["GeneSymbol"], pos, wt, mut)
            if sig in PATHOGENIC:
                path_keys.add(key)
            else:
                benign_keys.add(key)
    print(f"ClinVar pathogenic keys: {len(path_keys)}, benign keys: {len(benign_keys)}")
    return path_keys, benign_keys


def load_cosmic_keys(rbp_genes):
    keys = set()
    n_missense = 0
    with open(CMC_CHANGES, newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if row["category"] != "missense":
                continue
            n_missense += 1
            gene = row["GENE_NAME"]
            if gene not in rbp_genes:
                continue
            try:
                d = json.loads(row["details"])
                wt, mut, pos = d["wt"], d["mut"], d["pos"]
            except (json.JSONDecodeError, KeyError):
                continue
            if wt not in "ACDEFGHIKLMNPQRSTVWY" or mut not in "ACDEFGHIKLMNPQRSTVWY":
                continue
            keys.add((gene, pos, wt, mut))
    print(f"Total COSMIC missense variants: {n_missense}")
    print(f"COSMIC missense keys matching an RBP gene: {len(keys)}")
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
    protgps_wt, deeploc_wt = {}, {}
    seen = set()
    with open(MERGED_TABLE, newline="") as f:
        for row in csv.DictReader(f):
            acc = row["uniprot_accession"]
            if acc in seen:
                continue
            seen.add(acc)
            pvals, ok = [], True
            for c in PROTGPS_COLS:
                v = row.get(f"protgps_{c}", "").strip()
                if not v or v == "NA":
                    ok = False
                    break
                pvals.append(float(v))
            if ok:
                protgps_wt[acc] = pvals
            dvals, ok = [], True
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
    rbp_genes = load_rbp_genes()
    clinvar_path, clinvar_benign = load_clinvar_keys()
    cosmic_keys = load_cosmic_keys(rbp_genes)
    wt_tango = load_wt_tango()
    wt_protgps, wt_deeploc = load_wt_protgps_deeploc()
    print(f"WT baselines: {len(wt_tango)} TANGO, {len(wt_protgps)} protGPS, {len(wt_deeploc)} DeepLoc")

    buckets = {b: {"tango": [], "protgps": [], "deeploc": []}
               for b in ["clinvar_pathogenic", "clinvar_benign", "cosmic", "background"]}

    SAMPLE_STRIDE = 60
    n_total = 0
    counts_matched = {"clinvar_pathogenic": 0, "clinvar_benign": 0, "cosmic": 0}

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

            row_buckets = []
            if key in clinvar_path:
                row_buckets.append("clinvar_pathogenic")
            if key in clinvar_benign:
                row_buckets.append("clinvar_benign")
            if key in cosmic_keys:
                row_buckets.append("cosmic")
            is_background_sample = (i % SAMPLE_STRIDE == 0)
            if not row_buckets and not is_background_sample:
                continue

            tango_d = protgps_d = deeploc_d = None
            wt_agg = wt_tango.get(acc)
            if wt_agg is not None:
                try:
                    tango_d = abs(float(row["AGG"]) - wt_agg)
                except ValueError:
                    pass
            wt_p = wt_protgps.get(acc)
            if wt_p is not None:
                try:
                    protgps_d = sum(abs(float(row[c]) - w) for c, w in zip(PROTGPS_COLS, wt_p))
                except ValueError:
                    pass
            wt_d = wt_deeploc.get(acc)
            if wt_d is not None:
                try:
                    deeploc_d = sum(abs(float(row[c]) - w) for c, w in zip(DEEPLOC_NUMERIC_COLS, wt_d))
                except ValueError:
                    pass

            for b in row_buckets:
                counts_matched[b] += 1
                if tango_d is not None:
                    buckets[b]["tango"].append(tango_d)
                if protgps_d is not None:
                    buckets[b]["protgps"].append(protgps_d)
                if deeploc_d is not None:
                    buckets[b]["deeploc"].append(deeploc_d)

            if is_background_sample:
                if tango_d is not None:
                    buckets["background"]["tango"].append(tango_d)
                if protgps_d is not None:
                    buckets["background"]["protgps"].append(protgps_d)
                if deeploc_d is not None:
                    buckets["background"]["deeploc"].append(deeploc_d)

    print(f"\nTotal satmut rows: {n_total}")
    print(f"Matched counts: {counts_matched}")
    for b in buckets:
        print(f"  {b}: tango n={len(buckets[b]['tango'])}, protgps n={len(buckets[b]['protgps'])}, deeploc n={len(buckets[b]['deeploc'])}")

    with open(OUT_JSON, "w") as f:
        json.dump({"n_total": n_total, "counts_matched": counts_matched, "buckets": buckets}, f)
    print(f"\nWrote {OUT_JSON}")


if __name__ == "__main__":
    main()
