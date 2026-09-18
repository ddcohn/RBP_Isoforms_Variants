import csv
import sys

csv.field_size_limit(sys.maxsize)

TARGET = "/u/project/kappel/ddcohn/RBP/table_260823_with_rna.csv"
SOURCE = "/u/project/kappel/RBP/Misc/Classical_RBDs_Annotated.csv"

NEW_COLUMNS = [
    "Domains", "Domains_count", "Domains_avg_size", "Domains_total_size",
    "Domains_range", "Domains_discrete_seq",
    "Has_RBD", "RBD_name_ranges", "RBD_all_ranges", "RBD_names",
]

print("Loading + deduping source (Classical_RBDs_Annotated.csv)...")
source_data = {}
n_source_rows = 0
with open(SOURCE, newline="", encoding="utf-8", errors="replace") as f:
    reader = csv.DictReader(f)
    for row in reader:
        n_source_rows += 1
        uid = row.get("uniprot_id", "")
        if uid and uid not in source_data:
            source_data[uid] = {c: row.get(c, "") for c in NEW_COLUMNS}
        if n_source_rows % 5000 == 0:
            print(f"  ... {n_source_rows} source rows read", flush=True)

print(f"  {n_source_rows} source rows -> {len(source_data)} distinct proteins")

print("Reading target table...")
rows = []
with open(TARGET, newline="", encoding="utf-8", errors="replace") as f:
    reader = csv.DictReader(f)
    fieldnames = reader.fieldnames
    for row in reader:
        rows.append(row)
print(f"  {len(rows)} rows, {len(fieldnames)} existing columns")

for col in NEW_COLUMNS:
    if col in fieldnames:
        print(f"  WARNING: column {col} already exists, will be overwritten")

out_fieldnames = fieldnames + [c for c in NEW_COLUMNS if c not in fieldnames]

n_matched = 0
for row in rows:
    uid = row["uniprot_accession"]
    src = source_data.get(uid)
    if src:
        n_matched += 1
        for c in NEW_COLUMNS:
            row[c] = src[c]
    else:
        for c in NEW_COLUMNS:
            row[c] = ""

print(f"\nRows matched to Domains/RBD data: {n_matched}/{len(rows)} ({100*n_matched/len(rows):.1f}%)")

print(f"\nWriting merged table to {TARGET} ...")
with open(TARGET, "w", newline="", encoding="utf-8") as f:
    writer = csv.DictWriter(f, fieldnames=out_fieldnames)
    writer.writeheader()
    writer.writerows(rows)
print("Done.")
