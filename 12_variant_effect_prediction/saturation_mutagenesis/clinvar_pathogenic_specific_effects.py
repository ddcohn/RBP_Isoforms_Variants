"""
For just the ClinVar Pathogenic/Likely pathogenic missense variants
that match a real satmut mutation, answers two concrete questions:
  1. How many actually INCREASE TANGO's aggregation score vs WT
     (not just |delta| > 0 -- a real directional increase)?
  2. How many actually change DeepLoc's predicted Localizations
     category vs WT (a categorical flip, not just a numeric wobble)?
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
OUT_JSON = f"{BASE}/clinvar_pathogenic_specific_effects.json"

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


def load_wt_tango():
    wt = {}
    with open(WT_TANGO, newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            try:
                wt[row["uniprot_accession"]] = float(row["AGG"])
            except ValueError:
                pass
    return wt


def load_wt_deeploc_loc():
    wt = {}
    seen = set()
    with open(MERGED_TABLE, newline="") as f:
        for row in csv.DictReader(f):
            acc = row["uniprot_accession"]
            if acc in seen:
                continue
            seen.add(acc)
            loc = row.get("deeploc_Localizations", "").strip()
            if loc:
                wt[acc] = loc
    return wt


def main():
    pathogenic_keys = load_pathogenic_keys()
    wt_tango = load_wt_tango()
    wt_loc = load_wt_deeploc_loc()

    n_matched = 0
    n_tango_scored = 0
    n_agg_increased = 0
    n_agg_decreased = 0
    n_agg_unchanged = 0
    n_deeploc_scored = 0
    n_loc_changed = 0
    examples_agg_increase = []
    examples_loc_changed = []

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

            wt_agg = wt_tango.get(acc)
            if wt_agg is not None:
                try:
                    mut_agg = float(row["AGG"])
                    n_tango_scored += 1
                    diff = mut_agg - wt_agg
                    if diff > 0:
                        n_agg_increased += 1
                        if len(examples_agg_increase) < 10:
                            examples_agg_increase.append([row["mutation_id"], wt_agg, mut_agg, diff])
                    elif diff < 0:
                        n_agg_decreased += 1
                    else:
                        n_agg_unchanged += 1
                except ValueError:
                    pass

            wt_l = wt_loc.get(acc)
            if wt_l is not None:
                n_deeploc_scored += 1
                mut_l = row["Localizations"].strip()
                if mut_l != wt_l:
                    n_loc_changed += 1
                    if len(examples_loc_changed) < 10:
                        examples_loc_changed.append([row["mutation_id"], wt_l, mut_l])

    print(f"\nTotal satmut rows matching a Pathogenic/Likely pathogenic ClinVar missense variant: {n_matched}")
    print(f"\n--- TANGO aggregation ---")
    print(f"Scored (had WT baseline): {n_tango_scored}")
    print(f"Increased aggregation (mutant AGG > WT AGG): {n_agg_increased} ({100*n_agg_increased/n_tango_scored:.1f}%)")
    print(f"Decreased aggregation: {n_agg_decreased} ({100*n_agg_decreased/n_tango_scored:.1f}%)")
    print(f"Unchanged: {n_agg_unchanged}")
    print(f"\n--- DeepLoc localization ---")
    print(f"Scored (had WT baseline): {n_deeploc_scored}")
    print(f"Changed predicted localization category vs WT: {n_loc_changed} ({100*n_loc_changed/n_deeploc_scored:.1f}%)")

    out = {
        "n_matched": n_matched,
        "tango": {"n_scored": n_tango_scored, "n_increased": n_agg_increased,
                  "n_decreased": n_agg_decreased, "n_unchanged": n_agg_unchanged,
                  "examples_increase": examples_agg_increase},
        "deeploc": {"n_scored": n_deeploc_scored, "n_changed": n_loc_changed,
                    "examples_changed": examples_loc_changed},
    }
    with open(OUT_JSON, "w") as f:
        json.dump(out, f, indent=2)
    print(f"\nWrote {OUT_JSON}")


if __name__ == "__main__":
    main()
