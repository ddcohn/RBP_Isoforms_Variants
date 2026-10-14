import csv

D = "/Users/danielcohn/Desktop/RBP_repo_work/RBP_seq_STRING_OpenTargets/12_variant_effect_prediction/rbp_filter"
TABLE = f"{D}/rbp_wt_list_with_go.tsv"

rows_out = []
n_trna = n_rrna = 0
with open(TABLE, newline="") as f:
    r = csv.DictReader(f, delimiter="\t")
    fieldnames = r.fieldnames + ["tRNA_binding", "rRNA_binding"]
    for row in r:
        mf = row.get("GO_Molecular_Function", "")
        is_trna = "GO:0000049" in mf
        is_rrna = "GO:0019843" in mf
        row["tRNA_binding"] = "Y" if is_trna else "N"
        row["rRNA_binding"] = "Y" if is_rrna else "N"
        if is_trna:
            n_trna += 1
        if is_rrna:
            n_rrna += 1
        rows_out.append(row)

print(f"Total genes: {len(rows_out)}")
print(f"tRNA binding (GO:0000049): {n_trna}")
print(f"rRNA binding (GO:0019843): {n_rrna}")

with open(TABLE, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
    w.writeheader()
    w.writerows(rows_out)
print(f"\nUpdated {TABLE} in place")
