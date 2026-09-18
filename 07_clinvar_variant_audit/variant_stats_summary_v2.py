import csv
import sys
from collections import defaultdict

csv.field_size_limit(sys.maxsize)

TARGET = "/u/project/kappel/tchhabri/clinvar/Isoform_Post_Merge_PSLab_OpenTargets_withVariantStats.csv"

COLS = [
    "Benign_Not_ClassicalRBD", "Benign_In_ClassicalRBD",
    "Pathogenic_Not_ClassicalRBD", "Pathogenic_In_ClassicalRBD",
    "VUS_Not_ClassicalRBD", "VUS_In_ClassicalRBD",
    "Total_Pathogenic", "Total_VUS", "Total_Benign",
]


def to_num(v):
    v = (v or "").strip()
    if not v:
        return None
    try:
        return float(v)
    except ValueError:
        return None


# pass 1: collect ALL rows per protein, and check whether stat columns are
# actually consistent across every duplicate row for that protein (they
# should be, if these are meant to be per-protein rollups broadcast onto
# isoform rows)
per_protein_rows = defaultdict(list)
n_rows = 0
with open(TARGET, newline="", encoding="utf-8", errors="replace") as f:
    reader = csv.DictReader(f)
    for row in reader:
        n_rows += 1
        uid = row["uniprot_id"]
        vals = tuple(row.get(c) for c in COLS)
        per_protein_rows[uid].append(vals)
        if n_rows % 10000 == 0:
            print(f"  ... {n_rows} rows read", flush=True)

n_proteins = len(per_protein_rows)
print(f"\nTotal rows: {n_rows}")
print(f"Distinct proteins (uniprot_id): {n_proteins}")
print(f"Average rows per protein: {n_rows/n_proteins:.2f}")

dup_counts = sorted(((len(v), k) for k, v in per_protein_rows.items()), reverse=True)
print(f"\nTop 10 most-duplicated proteins:")
for count, uid in dup_counts[:10]:
    print(f"  {uid}: {count} rows")

# check consistency: for each protein, are all duplicate rows' stat columns identical?
n_inconsistent = 0
example_inconsistent = None
for uid, rows in per_protein_rows.items():
    if len(set(rows)) > 1:
        n_inconsistent += 1
        if example_inconsistent is None:
            example_inconsistent = (uid, rows[:3])

print(f"\nProteins where duplicate rows DISAGREE on stat columns: {n_inconsistent}/{n_proteins}")
if example_inconsistent:
    print(f"  Example: {example_inconsistent}")

# pass 2: dedupe (one row per protein, first occurrence) and compute real stats
sums = {c: 0.0 for c in COLS}
n_with_pathogenic = n_with_vus = n_with_benign = n_with_only_vus = n_with_any = 0
for uid, rows in per_protein_rows.items():
    vals = rows[0]  # first occurrence
    parsed = {c: (to_num(v) or 0.0) for c, v in zip(COLS, vals)}
    for c in COLS:
        sums[c] += parsed[c]
    tp, tv, tb = parsed["Total_Pathogenic"], parsed["Total_VUS"], parsed["Total_Benign"]
    if tp + tv + tb > 0:
        n_with_any += 1
    if tp > 0:
        n_with_pathogenic += 1
    if tv > 0:
        n_with_vus += 1
    if tb > 0:
        n_with_benign += 1
    if tv > 0 and tp == 0 and tb == 0:
        n_with_only_vus += 1

print(f"\n=== DEDUPLICATED (one row per protein) ===")
print(f"Column sums across {n_proteins} distinct proteins:")
for c in COLS:
    print(f"  {c}: {sums[c]:,.0f}")

total_variants = sums["Total_Pathogenic"] + sums["Total_VUS"] + sums["Total_Benign"]
print(f"\nOverall breakdown (of {total_variants:,.0f} total classified variant instances):")
print(f"  Pathogenic: {sums['Total_Pathogenic']:,.0f} ({100*sums['Total_Pathogenic']/total_variants:.1f}%)")
print(f"  VUS:        {sums['Total_VUS']:,.0f} ({100*sums['Total_VUS']/total_variants:.1f}%)")
print(f"  Benign:     {sums['Total_Benign']:,.0f} ({100*sums['Total_Benign']/total_variants:.1f}%)")

print(f"\nProtein-level coverage (out of {n_proteins} distinct proteins):")
print(f"  With >=1 classified variant of any kind: {n_with_any} ({100*n_with_any/n_proteins:.1f}%)")
print(f"  With >=1 Pathogenic: {n_with_pathogenic} ({100*n_with_pathogenic/n_proteins:.1f}%)")
print(f"  With >=1 VUS: {n_with_vus} ({100*n_with_vus/n_proteins:.1f}%)")
print(f"  With >=1 Benign: {n_with_benign} ({100*n_with_benign/n_proteins:.1f}%)")
print(f"  ONLY VUS (no path/benign at all): {n_with_only_vus} ({100*n_with_only_vus/n_proteins:.1f}%)")
