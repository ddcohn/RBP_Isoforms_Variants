"""
Same question as clinvar_pathogenic_specific_effects.py but for protGPS:
for ClinVar Pathogenic/Likely pathogenic missense variants, how many
change the protein's PREDICTED DOMINANT condensate (the single
compartment with the highest score) vs WT -- the categorical analog
of DeepLoc's "localization category changed" question.
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
MERGED_TABLE = "/u/project/kappel/ddcohn/RBP/merged_260903B_ptm_pred_slim_RBPs_domains_spurs_rsa_iso_hpa_with_rnp_rna_class.csv"
OUT_JSON = f"{BASE}/clinvar_pathogenic_protgps_effects.json"

PROTGPS_COLS = [
    "NUCLEAR_SPECKLE_Score", "P-BODY_Score", "PML-BDOY_Score",
    "POST_SYNAPTIC_DENSITY_Score", "STRESS_GRANULE_Score", "CHROMOSOME_Score",
    "NUCLEOLUS_Score", "NUCLEAR_PORE_COMPLEX_Score", "CAJAL_BODY_Score",
    "RNA_GRANULE_Score", "CELL_JUNCTION_Score", "TRANSCRIPTIONAL_Score",
]

AA3_TO_1 = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q",
    "Glu": "E", "Gly": "G", "His": "H", "Ile": "I", "Leu": "L", "Lys": "K",
    "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W",
    "Tyr": "Y", "Val": "V",
}
MUT_ID_RE = re.compile(r"^(?P<acc>[^|]+)\|(?P<gene>[^|]+)\|(?P<wt>[A-Z])(?P<pos>\d+)(?P<mut>[A-Z])$")
PATHOGENIC = {"Pathogenic", "Likely pathogenic", "Pathogenic/Likely pathogenic"}


def load_pathogenic_keys():
    var_sig = {}
    with open(CLINVAR_META_RBP, newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            var_sig[row["VariationID"]] = row["ClinicalSignificance"]
    keys = set()
    with open(CLINVAR_CHANGES, newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if row["category"] != "missense":
                continue
            sig = var_sig.get(row["VariationID"])
            if sig not in PATHOGENIC:
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
            keys.add((row["GeneSymbol"], pos, wt, mut))
    print(f"Pathogenic/Likely pathogenic missense keys: {len(keys)}")
    return keys


def load_wt_protgps():
    wt = {}
    seen = set()
    with open(MERGED_TABLE, newline="") as f:
        for row in csv.DictReader(f):
            acc = row["uniprot_accession"]
            if acc in seen:
                continue
            seen.add(acc)
            vals = []
            ok = True
            for c in PROTGPS_COLS:
                v = row.get(f"protgps_{c}", "").strip()
                if not v or v == "NA":
                    ok = False
                    break
                vals.append(float(v))
            if ok:
                wt[acc] = vals
    return wt


def dominant(vals):
    best_i = max(range(len(vals)), key=lambda i: vals[i])
    return PROTGPS_COLS[best_i], vals[best_i]


def main():
    pathogenic_keys = load_pathogenic_keys()
    wt_protgps = load_wt_protgps()

    n_matched = 0
    n_scored = 0
    n_dominant_changed = 0
    examples_changed = []

    with open(SATMUT_ALL, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            m = MUT_ID_RE.match(row["mutation_id"])
            if not m:
                continue
            acc = m.group("acc")
            gene = m.group("gene")
            pos = int(m.group("pos"))
            wt_aa = m.group("wt")
            mut_aa = m.group("mut")
            key = (gene, pos, wt_aa, mut_aa)
            if key not in pathogenic_keys:
                continue
            n_matched += 1

            wt_vals = wt_protgps.get(acc)
            if wt_vals is None:
                continue
            try:
                mut_vals = [float(row[c]) for c in PROTGPS_COLS]
            except ValueError:
                continue
            n_scored += 1

            wt_dom, wt_dom_score = dominant(wt_vals)
            mut_dom, mut_dom_score = dominant(mut_vals)
            if wt_dom != mut_dom:
                n_dominant_changed += 1
                if len(examples_changed) < 10:
                    examples_changed.append([
                        row["mutation_id"], wt_dom, round(wt_dom_score, 4),
                        mut_dom, round(mut_dom_score, 4),
                    ])

    print(f"\nTotal satmut rows matching a Pathogenic/Likely pathogenic ClinVar missense variant: {n_matched}")
    print(f"\n--- protGPS dominant condensate ---")
    print(f"Scored (had WT baseline): {n_scored}")
    print(f"Changed predicted dominant condensate vs WT: {n_dominant_changed} ({100*n_dominant_changed/n_scored:.1f}%)")

    out = {"n_matched": n_matched, "n_scored": n_scored, "n_dominant_changed": n_dominant_changed,
           "examples_changed": examples_changed}
    with open(OUT_JSON, "w") as f:
        json.dump(out, f, indent=2)
    print(f"\nWrote {OUT_JSON}")


if __name__ == "__main__":
    main()
