"""
Adds ENCORI's top-enriched-motif data onto the merged table (after
add_encori_clip.py), joined by gene_symbol only -- not a compound
join against the per-cell-line rows already present, because a
motif's own `cellTissue` field is itself a pipe-separated list of
multiple cell lines (e.g. "BE(2)-C|HEK293|HeLa|HepG2|K562"), not a
single value that could align 1:1 with one binding-site row. This is
a genuine cross-product: a gene with N existing rows (one per binding-
site cell line) and M motif rows becomes N*M rows. Motifs are
gene-level findings (the top motif applies across all of a gene's
experiments), so this isn't incoherent, just real further row growth
-- confirmed acceptable before running.
"""
import csv

csv.field_size_limit(2**31 - 1)

TABLE = "/u/project/kappel/ddcohn/RBP/merged_260903B_ptm_pred_slim_RBPs_domains_spurs_rsa_iso_hpa_with_rnp_rna_class.csv"
MOTIFS = "/u/project/kappel/ddcohn/protein_variant_effects/rbp_only/encori_motifs_rbp_filtered.tsv"

with open(MOTIFS, newline="") as f:
    r = csv.DictReader(f, delimiter="\t")
    motif_cols = [c for c in r.fieldnames if c != "RBP"]
    motifs_by_gene = {}
    for row in r:
        motifs_by_gene.setdefault(row["RBP"], []).append(
            {f"encori_motif_{c}": row[c] for c in motif_cols}
        )

print(f"Motifs: {len(motifs_by_gene)} genes, {sum(len(v) for v in motifs_by_gene.values())} rows loaded")

with open(TABLE, newline="") as f:
    r = csv.DictReader(f)
    orig_fieldnames = r.fieldnames
    new_cols = [f"encori_motif_{c}" for c in motif_cols]
    assert not any(c in orig_fieldnames for c in new_cols), "unexpected column collision"
    fieldnames = orig_fieldnames + new_cols
    rows = list(r)

out_rows = []
n_matched_rows = 0
n_extra_rows = 0
for row in rows:
    gene = row["gene_symbol"]
    matches = motifs_by_gene.get(gene)
    if not matches:
        new_row = dict(row)
        for c in new_cols:
            new_row[c] = ""
        out_rows.append(new_row)
    else:
        n_matched_rows += 1
        n_extra_rows += len(matches) - 1
        for m in matches:
            new_row = dict(row)
            new_row.update(m)
            out_rows.append(new_row)

print(f"Table rows in: {len(rows)}")
print(f"Rows matched to >=1 motif: {n_matched_rows}")
print(f"Extra rows from cross-product (>1 motif per gene): {n_extra_rows}")
print(f"Table rows out: {len(out_rows)}")

with open(TABLE, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=fieldnames)
    w.writeheader()
    w.writerows(out_rows)

print(f"\nUpdated {TABLE} in place")
