"""
Adds this project's WT-baseline protGPS and DeepLoc predictions onto
the asharma-table-plus-RNP/RNA-class merge (see merge_into_asharma_table.py),
joined by UniProt accession, column names prefixed protgps_/deeploc_ to
keep sources unambiguous. One row per protein in, one row per protein
out -- no cell-line/multi-row complication like the ENCORI CLIP data has.

TANGO has no equivalent WT-baseline run (only ever run on the 31.5M
saturation mutants, not WT sequences) -- deliberately not included here,
can be added later if that baseline gets run.
"""
import csv

csv.field_size_limit(2**31 - 1)

TABLE = "/u/project/kappel/ddcohn/RBP/merged_260903B_ptm_pred_slim_RBPs_domains_spurs_rsa_iso_hpa_with_rnp_rna_class.csv"
PROTGPS = "/u/project/kappel/ddcohn/protein_variant_effects/rbp_only/rbp_wt_protgps_results.tsv"
DEEPLOC = "/u/project/kappel/ddcohn/protein_variant_effects/rbp_only/rbp_wt_deeploc_results.csv"

def load(path, id_col, delim, prefix):
    data = {}
    with open(path, newline="") as f:
        r = csv.DictReader(f, delimiter=delim)
        cols = [c for c in r.fieldnames if c != id_col]
        for row in r:
            acc = row[id_col].split("|", 1)[1] if "|" in row[id_col] else row[id_col]
            data[acc] = {f"{prefix}{c}": row[c] for c in cols}
    return cols, data

pg_cols, pg_data = load(PROTGPS, "id", "\t", "protgps_")
dl_cols, dl_data = load(DEEPLOC, "Protein_ID", ",", "deeploc_")

print(f"protGPS WT: {len(pg_data)} proteins loaded")
print(f"DeepLoc WT: {len(dl_data)} proteins loaded")

with open(TABLE, newline="") as f:
    r = csv.DictReader(f)
    orig_fieldnames = r.fieldnames
    new_cols = [f"protgps_{c}" for c in pg_cols] + [f"deeploc_{c}" for c in dl_cols]
    assert not any(c in orig_fieldnames for c in new_cols), "unexpected column collision"
    fieldnames = orig_fieldnames + new_cols
    rows = list(r)

n_pg_matched = n_dl_matched = 0
for row in rows:
    acc = row["uniprot_accession"]
    pg = pg_data.get(acc)
    dl = dl_data.get(acc)
    if pg:
        row.update(pg)
        n_pg_matched += 1
    else:
        for c in pg_cols:
            row[f"protgps_{c}"] = ""
    if dl:
        row.update(dl)
        n_dl_matched += 1
    else:
        for c in dl_cols:
            row[f"deeploc_{c}"] = ""

print(f"Table rows: {len(rows)}")
print(f"Matched protGPS: {n_pg_matched}")
print(f"Matched DeepLoc: {n_dl_matched}")

with open(TABLE, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=fieldnames)
    w.writeheader()
    w.writerows(rows)

print(f"\nUpdated {TABLE} in place")
