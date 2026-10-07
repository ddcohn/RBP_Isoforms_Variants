"""
Extracts the "Binding Site Summary" sheet from the user's ENCORI reference
workbook (CLIP-seq-derived, per Zhou et al. 2026 Nature Methods) as-is --
one row per RBP x cell line/tissue, no averaging or collapsing across
experiments. Filtered to the RBPs in this project's 2,335-gene list.

A separate table, not joined onto the one-row-per-protein isoform table:
79 of the 209 ENCORI-covered RBPs have multiple cell-line rows, and
averaging them would discard real experiment-specific binding-region
data the user explicitly wants kept.
"""
import csv
import openpyxl

XLSX = "/Users/danielcohn/Desktop/encori/ENCORI_human_RBP_reference.xlsx"
RBP_LIST = "/Users/danielcohn/Desktop/RBP_repo_work/RBP_seq_STRING_OpenTargets/12_variant_effect_prediction/rbp_filter/rbp_wt_list_with_go.tsv"
OUT = "/Users/danielcohn/Desktop/RBP_repo_work/RBP_seq_STRING_OpenTargets/12_variant_effect_prediction/rbp_filter/encori_binding_site_summary_rbp_filtered.tsv"

rbp_symbols = set()
with open(RBP_LIST, newline="") as f:
    r = csv.DictReader(f, delimiter="\t")
    for row in r:
        rbp_symbols.add(row["GeneSymbol"])

print(f"RBP list: {len(rbp_symbols)} gene symbols loaded")

wb = openpyxl.load_workbook(XLSX, read_only=True)
ws = wb["Binding Site Summary"]
rows = list(ws.iter_rows(min_row=1, values_only=True))
header = rows[2]
data = rows[3:]

matched = [row for row in data if row[0] in rbp_symbols]
unique_rbps = set(row[0] for row in matched)

print(f"ENCORI sheet total rows: {len(data)}")
print(f"ENCORI unique RBPs: {len(set(row[0] for row in data))}")
print(f"Matched rows (RBP in our list): {len(matched)}")
print(f"Matched unique RBPs: {len(unique_rbps)}")

with open(OUT, "w", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(header)
    w.writerows(matched)

print(f"Wrote {OUT}")
