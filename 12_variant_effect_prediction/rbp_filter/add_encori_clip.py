"""
Adds ENCORI CLIP-seq binding-region data onto the merged table (see
merge_into_asharma_table.py / add_protgps_deeploc_wt.py), joined by
gene_symbol. Unlike those two scripts, this one is NOT a clean 1:1
join -- 79 of the 195 ENCORI-matched genes have more than one row
(one per cell line/tissue profiled), and per explicit instruction not
to average across cell lines, those genes' entire row (all 187
preceding columns) gets duplicated once per ENCORI row rather than
collapsed. Genes with zero or one ENCORI row are unaffected and stay
single-row. Table grows from 2,332 to 2,459 rows as a result -- expect
duplicates on gene_symbol/uniprot_accession for any multi-cell-line gene.
"""
import csv

csv.field_size_limit(2**31 - 1)

TABLE = "/u/project/kappel/ddcohn/RBP/merged_260903B_ptm_pred_slim_RBPs_domains_spurs_rsa_iso_hpa_with_rnp_rna_class.csv"
ENCORI = "/u/project/kappel/ddcohn/protein_variant_effects/rbp_only/encori_binding_site_summary_rbp_filtered.tsv"

with open(ENCORI, newline="") as f:
    r = csv.DictReader(f, delimiter="\t")
    encori_cols = [c for c in r.fieldnames if c != "RBP"]
    encori_by_gene = {}
    for row in r:
        encori_by_gene.setdefault(row["RBP"], []).append(
            {f"encori_{c}": row[c] for c in encori_cols}
        )

print(f"ENCORI: {len(encori_by_gene)} genes, {sum(len(v) for v in encori_by_gene.values())} rows loaded")

with open(TABLE, newline="") as f:
    r = csv.DictReader(f)
    orig_fieldnames = r.fieldnames
    new_cols = [f"encori_{c}" for c in encori_cols]
    assert not any(c in orig_fieldnames for c in new_cols), "unexpected column collision"
    fieldnames = orig_fieldnames + new_cols
    rows = list(r)

out_rows = []
n_matched_genes = 0
n_extra_rows = 0
for row in rows:
    gene = row["gene_symbol"]
    matches = encori_by_gene.get(gene)
    if not matches:
        new_row = dict(row)
        for c in new_cols:
            new_row[c] = ""
        out_rows.append(new_row)
    else:
        n_matched_genes += 1
        n_extra_rows += len(matches) - 1
        for m in matches:
            new_row = dict(row)
            new_row.update(m)
            out_rows.append(new_row)

print(f"Table rows in: {len(rows)}")
print(f"Genes matched to ENCORI: {n_matched_genes}")
print(f"Extra duplicate rows from multi-cell-line genes: {n_extra_rows}")
print(f"Table rows out: {len(out_rows)}")

with open(TABLE, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=fieldnames)
    w.writeheader()
    w.writerows(out_rows)

print(f"\nUpdated {TABLE} in place")
