import csv
import sys

csv.field_size_limit(sys.maxsize)

OLD = "/u/project/kappel/ddcohn/RBP/table_260823_with_rna.csv.bak_pre_domains_rbd"
NEW = "/u/project/kappel/ddcohn/RBP/table_260823_with_rna.csv"


def load(path):
    with open(path, newline="", encoding="utf-8", errors="replace") as f:
        reader = csv.reader(f)
        header = next(reader)
        rows = list(reader)
    return header, rows


old_header, old_rows = load(OLD)
new_header, new_rows = load(NEW)

print(f"Old: {len(old_rows)} rows, {len(old_header)} columns")
print(f"New: {len(new_rows)} rows, {len(new_header)} columns")

new_cols_added = [c for c in new_header if c not in old_header]
print(f"New columns added: {new_cols_added}")
missing_old_cols = [c for c in old_header if c not in new_header]
print(f"Old columns missing in new (should be empty): {missing_old_cols}")

if len(old_rows) != len(new_rows):
    print(f"MISMATCH: row count differs! old={len(old_rows)} new={len(new_rows)}")
else:
    print("Row count matches.")

ragged_old = sum(1 for r in old_rows if len(r) != len(old_header))
ragged_new = sum(1 for r in new_rows if len(r) != len(new_header))
print(f"Ragged rows -- old: {ragged_old}, new: {ragged_new}")

old_uid_idx = old_header.index("uniprot_accession")
new_uid_idx = new_header.index("uniprot_accession")
old_by_uid = {r[old_uid_idx]: r for r in old_rows}
new_by_uid = {r[new_uid_idx]: r for r in new_rows}

missing_in_new = set(old_by_uid) - set(new_by_uid)
extra_in_new = set(new_by_uid) - set(old_by_uid)
print(f"UIDs in old but missing from new: {len(missing_in_new)}")
print(f"UIDs in new but not in old: {len(extra_in_new)}")

shared_cols = [c for c in old_header if c in new_header]
mismatches = {c: 0 for c in shared_cols}
sample_mismatches = {}
checked = 0
for uid, old_row in old_by_uid.items():
    new_row = new_by_uid.get(uid)
    if not new_row:
        continue
    checked += 1
    for c in shared_cols:
        oi = old_header.index(c)
        ni = new_header.index(c)
        if old_row[oi] != new_row[ni]:
            mismatches[c] += 1
            if c not in sample_mismatches:
                sample_mismatches[c] = (uid, old_row[oi][:80], new_row[ni][:80])

print(f"\nChecked {checked} rows for old-column integrity.")
any_mismatch = False
for c, n in mismatches.items():
    if n > 0:
        any_mismatch = True
        print(f"  MISMATCH in column '{c}': {n} rows differ. Example: {sample_mismatches[c]}")
if not any_mismatch:
    print("All shared (pre-existing) columns are byte-identical between old and new, for every row. No corruption/misalignment.")

# quick coverage summary of new columns
n_has_domains = sum(1 for r in new_rows if r[new_header.index("Domains")].strip() not in ("", "[]"))
n_has_rbd = sum(1 for r in new_rows if r[new_header.index("Has_RBD")].strip() == "True" or r[new_header.index("Has_RBD")].strip() == "1")
print(f"\nRows with non-empty Domains: {n_has_domains}/{len(new_rows)}")
print(f"Rows with Has_RBD flagged true: {n_has_rbd}/{len(new_rows)}")
