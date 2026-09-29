import csv
import sys

csv.field_size_limit(sys.maxsize)
D = "/u/project/kappel/ddcohn/protein_variant_effects/rbp_only"

# uniprot_accession -> (GO_CC, GO_BP, GO_MF) from the user's own gene-level table
go = {}
n_rows = 0
with open("/u/project/kappel/ddcohn/RBP/table_260823_with_rna.csv", newline="", encoding="utf-8", errors="replace") as f:
    r = csv.DictReader(f)
    for row in r:
        n_rows += 1
        acc = row.get("uniprot_accession", "").strip()
        if acc:
            go[acc] = (
                row.get("GO_Cellular_Component", ""),
                row.get("GO_Biological_Process", ""),
                row.get("GO_Molecular_Function", ""),
            )
print(f"Rows scanned in table_260823_with_rna.csv: {n_rows}")
print(f"Unique UniProt accessions with GO data: {len(go)}")

n_total = 0
n_matched = 0
n_matched_nonempty = 0
rows_out = []
with open(f"{D}/rbp_wt_list.tsv", newline="") as f:
    r = csv.DictReader(f, delimiter="\t")
    fieldnames = r.fieldnames + ["GO_Cellular_Component", "GO_Biological_Process", "GO_Molecular_Function"]
    for row in r:
        n_total += 1
        acc = row["UniProtID"]
        cc, bp, mf = go.get(acc, ("", "", ""))
        if acc in go:
            n_matched += 1
            if cc or bp or mf:
                n_matched_nonempty += 1
        row["GO_Cellular_Component"] = cc
        row["GO_Biological_Process"] = bp
        row["GO_Molecular_Function"] = mf
        rows_out.append(row)

print(f"\nRBP genes: {n_total}")
print(f"Matched by UniProt accession: {n_matched}")
print(f"Matched with at least one non-empty GO field: {n_matched_nonempty}")
print(f"Unmatched: {n_total - n_matched}")

with open(f"{D}/rbp_wt_list_with_go.tsv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
    w.writeheader()
    w.writerows(rows_out)
print(f"\nWrote {D}/rbp_wt_list_with_go.tsv")
