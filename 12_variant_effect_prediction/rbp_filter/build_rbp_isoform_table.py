"""
Filters the full isoform table (table_260823_with_rna.csv, 20,431 proteins,
not RBP-specific) down to the 2,335 RBPs used throughout this project, and
joins in the HGNC_Gene_Family column already fetched for that RBP set.

A separate table from rbp_wt_list_with_go.tsv (which has a handful of
curated columns for the variant-effect pipeline specifically) -- this one
keeps every column of the real isoform table, just row-filtered + one
column added, so it's usable for isoform-level analyses that need more
than what the variant-effect gene list carries.
"""
import csv

csv.field_size_limit(2**31 - 1)

ISOFORM_TABLE = "/u/project/kappel/ddcohn/RBP/table_260823_with_rna.csv"
RBP_LIST = "/tmp/ddcohn_tmp/rbp_wt_list_with_go.tsv"
OUT = "/u/project/kappel/ddcohn/RBP/isoform_table_rbp_only_with_hgnc_family.csv"

family_by_uniprot = {}
with open(RBP_LIST, newline="") as f:
    r = csv.DictReader(f, delimiter="\t")
    for row in r:
        family_by_uniprot[row["UniProtID"]] = row["HGNC_Gene_Family"]

print(f"RBP list: {len(family_by_uniprot)} UniProt IDs loaded")

n_total = 0
n_matched = 0
rows_out = []
with open(ISOFORM_TABLE, newline="", encoding="utf-8", errors="replace") as f:
    r = csv.DictReader(f)
    fieldnames = r.fieldnames + ["HGNC_Gene_Family"]
    for row in r:
        n_total += 1
        acc = row["uniprot_accession"]
        if acc in family_by_uniprot:
            row["HGNC_Gene_Family"] = family_by_uniprot[acc]
            rows_out.append(row)
            n_matched += 1

print(f"Isoform table total rows scanned: {n_total}")
print(f"Matched to RBP list: {n_matched}")
print(f"Not matched (in RBP list but not in isoform table): {len(family_by_uniprot) - n_matched}")

with open(OUT, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=fieldnames)
    w.writeheader()
    w.writerows(rows_out)

print(f"Wrote {OUT}")
