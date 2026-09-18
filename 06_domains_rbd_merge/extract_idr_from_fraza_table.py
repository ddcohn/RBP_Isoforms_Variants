import csv
import sys
from collections import defaultdict

csv.field_size_limit(sys.maxsize)

PATH = "/u/project/kappel/RBP/Isoform_Table/Isoform_Post_Merge_PSLab_OpenTargets_Updated_Interim_20260817.csv"
OUT = "/u/home/d/ddcohn/_claude_idr_from_fraza.csv"

IDR_COLS = [
    "idr_method",
    "IDR_count", "IDR_avg_size", "IDR_total_size", "IDR_range",
    "IDR_discrete_seq", "IDR_concat_seq",
    "IDR_FCR", "IDR_NCPR", "IDR_isoelectric_point", "IDR_molecular_weight",
    "IDR_countNeg", "IDR_countPos", "IDR_countNeut",
    "IDR_fraction_negative", "IDR_fraction_positive", "IDR_fraction_expanding",
    "IDR_amino_acid_fractions", "IDR_fraction_disorder_promoting",
    "IDR_kappa", "IDR_mean_net_charge", "IDR_mean_hydropathy",
    "IDR_uversky_hydropathy", "IDR_PPII_propensity", "IDR_delta", "IDR_deltaMax",
]

per_protein_rows = defaultdict(list)
n_rows = 0
method_counts = defaultdict(int)

with open(PATH, newline="", encoding="utf-8", errors="replace") as f:
    reader = csv.DictReader(f)
    fieldnames = reader.fieldnames
    missing_cols = [c for c in IDR_COLS if c not in fieldnames]
    print(f"Columns not found in file: {missing_cols}")
    for row in reader:
        n_rows += 1
        uid = row.get("uniprot_id", "")
        vals = tuple(row.get(c, "") for c in IDR_COLS)
        per_protein_rows[uid].append(vals)
        method_counts[row.get("idr_method", "")] += 1
        if n_rows % 10000 == 0:
            print(f"  ... {n_rows} rows read", flush=True)

print(f"\nTotal rows: {n_rows}")
print(f"Distinct uniprot_id: {len(per_protein_rows)}")
print(f"\nidr_method value counts:")
for m, c in sorted(method_counts.items(), key=lambda x: -x[1]):
    print(f"  '{m}': {c}")

# consistency check across duplicate rows per protein
n_consistent = 0
n_inconsistent = 0
example_inconsistent = None
for uid, rows in per_protein_rows.items():
    if len(set(rows)) == 1:
        n_consistent += 1
    else:
        n_inconsistent += 1
        if example_inconsistent is None:
            example_inconsistent = (uid, rows[:3])

print(f"\nProteins with fully consistent IDR values across all their duplicate rows: {n_consistent}")
print(f"Proteins with INCONSISTENT IDR values across duplicate rows: {n_inconsistent}")
if example_inconsistent:
    print(f"Example inconsistent case: {example_inconsistent}")

# write deduplicated output (first occurrence per protein)
with open(OUT, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["uniprot_id"] + IDR_COLS)
    for uid, rows in per_protein_rows.items():
        writer.writerow([uid] + list(rows[0]))

print(f"\nWrote {len(per_protein_rows)} deduplicated rows to {OUT}")
