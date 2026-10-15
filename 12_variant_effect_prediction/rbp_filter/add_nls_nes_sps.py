"""
Adds the 32 NLS/NES/SPS columns from asharma's updated merged table
(merged_260903B_ptm_pred_slim_RBPs_domains_spurs_rsa_iso_hpa_nls_nes_sps.csv)
onto our already-merged table. Confirmed first that all other 137
columns in the source are identical to what's already merged in (see
merge_into_asharma_table.py) -- only these 32 are genuinely new:
classical NLS predictions (cNLS, mono- and bi-partite), NES
predictions, and a second independent method (SPS) for both NLS and
NES.

This is a per-protein broadcast join, not a row-multiplying one: our
table already has duplicate rows per protein from the ENCORI merges
(add_encori_clip.py / add_encori_motifs.py), and NLS/NES/SPS data is
protein-level (one set of values per uniprot_accession), so the same
values get repeated across all of a protein's existing duplicate rows
rather than creating new ones.
"""
import csv

csv.field_size_limit(2**31 - 1)

TABLE = "/u/project/kappel/ddcohn/RBP/merged_260903B_ptm_pred_slim_RBPs_domains_spurs_rsa_iso_hpa_with_rnp_rna_class.csv"
NEW_SOURCE = "/u/project/kappel/asharma/RBP/merged/merged_260903B_ptm_pred_slim_RBPs_domains_spurs_rsa_iso_hpa_nls_nes_sps.csv"

NEW_COLS = [
    "cNLS_mono_count", "cNLS_mono_coords", "cNLS_mono_scores", "cNLS_mono_seqs",
    "cNLS_bi_count", "cNLS_bi_coords", "cNLS_bi_scores", "cNLS_bi_seqs",
    "cNLS_count", "cNLS_regions", "cNLS_max_score",
    "NES_count", "NES_coords", "NES_scores", "NES_classes", "NES_seqs", "NES_max_score",
    "SPS_NLS_count", "SPS_NLS_coords", "SPS_NLS_SCR", "SPS_NLS_proba", "SPS_NLS_seqs",
    "SPS_NLS_candidates", "SPS_NLS_best_SCR",
    "SPS_NES_count", "SPS_NES_coords", "SPS_NES_SCR", "SPS_NES_proba", "SPS_NES_seqs",
    "SPS_NES_candidates", "SPS_NES_best_SCR", "SPS_status",
]

with open(NEW_SOURCE, newline="") as f:
    r = csv.DictReader(f)
    new_data = {row["uniprot_accession"]: {c: row[c] for c in NEW_COLS} for row in r}

print(f"NLS/NES/SPS source: {len(new_data)} proteins loaded")

with open(TABLE, newline="") as f:
    r = csv.DictReader(f)
    orig_fieldnames = r.fieldnames
    assert not any(c in orig_fieldnames for c in NEW_COLS), "unexpected column collision"
    fieldnames = orig_fieldnames + NEW_COLS
    rows = list(r)

n_matched = 0
for row in rows:
    acc = row["uniprot_accession"]
    extra = new_data.get(acc)
    if extra:
        row.update(extra)
        n_matched += 1
    else:
        for c in NEW_COLS:
            row[c] = ""

print(f"Table rows: {len(rows)}")
print(f"Matched to NLS/NES/SPS data: {n_matched}")
print(f"Not matched (left blank): {len(rows) - n_matched}")

with open(TABLE, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=fieldnames)
    w.writeheader()
    w.writerows(rows)

print(f"\nUpdated {TABLE} in place")
