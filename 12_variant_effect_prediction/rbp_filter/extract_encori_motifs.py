"""
Extracts the "Motifs" sheet from the user's ENCORI reference workbook,
filtered to this project's RBP list, same pattern and same no-averaging
rule as extract_encori_binding_summary.py: kept row-for-row as-is (up
to 2 rows per gene, "top 2 enriched motifs per RBP experiment" per the
sheet's own description -- not duplicates to collapse).
"""
import csv
import openpyxl

XLSX = "/Users/danielcohn/Desktop/encori/ENCORI_human_RBP_reference.xlsx"
RBP_LIST = "/Users/danielcohn/Desktop/RBP_repo_work/RBP_seq_STRING_OpenTargets/12_variant_effect_prediction/rbp_filter/rbp_wt_list_with_go.tsv"
OUT = "/Users/danielcohn/Desktop/RBP_repo_work/RBP_seq_STRING_OpenTargets/12_variant_effect_prediction/rbp_filter/encori_motifs_rbp_filtered.tsv"

rbp_symbols = set()
with open(RBP_LIST, newline="") as f:
    r = csv.DictReader(f, delimiter="\t")
    for row in r:
        rbp_symbols.add(row["GeneSymbol"])

print(f"RBP list: {len(rbp_symbols)} gene symbols loaded")

wb = openpyxl.load_workbook(XLSX, read_only=True)
ws = wb["Motifs"]
rows = list(ws.iter_rows(min_row=1, values_only=True))
header = rows[2]
data = rows[3:]

matched = [row for row in data if row[0] in rbp_symbols]
unique_rbps = set(row[0] for row in matched)

print(f"ENCORI Motifs sheet total rows: {len(data)}")
print(f"ENCORI Motifs unique RBPs: {len(set(row[0] for row in data))}")
print(f"Matched rows (RBP in our list): {len(matched)}")
print(f"Matched unique RBPs: {len(unique_rbps)}")

with open(OUT, "w", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(header)
    w.writerows(matched)

print(f"Wrote {OUT}")
