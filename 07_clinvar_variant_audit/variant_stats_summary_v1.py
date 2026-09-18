import csv
import sys

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
        return 0.0
    try:
        return float(v)
    except ValueError:
        return 0.0


n_rows = 0
sums = {c: 0.0 for c in COLS}
n_with_any_variant = 0
n_with_pathogenic = 0
n_with_vus = 0
n_with_benign = 0
n_with_only_vus = 0  # proteins whose ENTIRE variant record is VUS (no path/benign at all) -- most in need of resolution

with open(TARGET, newline="", encoding="utf-8", errors="replace") as f:
    reader = csv.DictReader(f)
    for row in reader:
        n_rows += 1
        vals = {}
        for c in COLS:
            v = to_num(row.get(c))
            sums[c] += v
            vals[c] = v
        tp, tv, tb = vals["Total_Pathogenic"], vals["Total_VUS"], vals["Total_Benign"]
        if tp + tv + tb > 0:
            n_with_any_variant += 1
        if tp > 0:
            n_with_pathogenic += 1
        if tv > 0:
            n_with_vus += 1
        if tb > 0:
            n_with_benign += 1
        if tv > 0 and tp == 0 and tb == 0:
            n_with_only_vus += 1
        if n_rows % 5000 == 0:
            print(f"  ... {n_rows} rows processed", flush=True)

print(f"\nTotal rows (proteins): {n_rows}")
print(f"\nRaw column sums (total variant-instances across all proteins):")
for c in COLS:
    print(f"  {c}: {sums[c]:,.0f}")

total_variants = sums["Total_Pathogenic"] + sums["Total_VUS"] + sums["Total_Benign"]
print(f"\nOverall breakdown (of {total_variants:,.0f} total classified variant instances):")
print(f"  Pathogenic: {sums['Total_Pathogenic']:,.0f} ({100*sums['Total_Pathogenic']/total_variants:.1f}%)")
print(f"  VUS:        {sums['Total_VUS']:,.0f} ({100*sums['Total_VUS']/total_variants:.1f}%)")
print(f"  Benign:     {sums['Total_Benign']:,.0f} ({100*sums['Total_Benign']/total_variants:.1f}%)")

print(f"\nProtein-level coverage (out of {n_rows} proteins):")
print(f"  Proteins with >=1 classified variant of any kind: {n_with_any_variant} ({100*n_with_any_variant/n_rows:.1f}%)")
print(f"  Proteins with >=1 Pathogenic variant: {n_with_pathogenic} ({100*n_with_pathogenic/n_rows:.1f}%)")
print(f"  Proteins with >=1 VUS: {n_with_vus} ({100*n_with_vus/n_rows:.1f}%)")
print(f"  Proteins with >=1 Benign variant: {n_with_benign} ({100*n_with_benign/n_rows:.1f}%)")
print(f"  Proteins whose ONLY classified variants are VUS (no path/benign at all): {n_with_only_vus} ({100*n_with_only_vus/n_rows:.1f}%)")

# RBD-split VUS breakdown specifically (relevant since pathogenicity-by-RBD is a stated analysis axis)
vus_in_rbd = sums["VUS_In_ClassicalRBD"]
vus_not_rbd = sums["VUS_Not_ClassicalRBD"]
print(f"\nVUS by classical RBD membership:")
print(f"  In classical RBD:  {vus_in_rbd:,.0f}")
print(f"  Not in classical RBD: {vus_not_rbd:,.0f}")
