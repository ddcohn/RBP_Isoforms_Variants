import csv
import json

D = "/u/project/kappel/ddcohn/protein_variant_effects/rbp_only"
csv.field_size_limit(2**31 - 1)

families = json.load(open(f"{D}/hgnc_families.json"))

rows_out = []
n_with_family = 0
with open(f"{D}/rbp_wt_list_with_go.tsv", newline="") as f:
    r = csv.DictReader(f, delimiter="\t")
    assert len(r.fieldnames) == 11, f"expected 11 columns, got {len(r.fieldnames)}: {r.fieldnames}"
    fieldnames = r.fieldnames + ["HGNC_Gene_Family"]
    for row in r:
        groups = families.get(row["GeneSymbol"]) or []
        row["HGNC_Gene_Family"] = ";".join(groups)
        if groups:
            n_with_family += 1
        rows_out.append(row)

print(f"Total genes: {len(rows_out)}")
print(f"With at least one HGNC family: {n_with_family}")
print(f"Output columns: {len(fieldnames)} -> {fieldnames}")

with open(f"{D}/rbp_wt_list_with_go.tsv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
    w.writeheader()
    w.writerows(rows_out)
print(f"\nUpdated {D}/rbp_wt_list_with_go.tsv in place")
