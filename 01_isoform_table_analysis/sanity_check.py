import csv
import sys
from collections import Counter, defaultdict

path = "/u/project/kappel/RBP/Isoform_Table/Isoform_Post_Merge_PSLab_OpenTargets_Updated_Interim_20260817.csv"

csv.field_size_limit(sys.maxsize)

total_rows = 0
uniprot_ids = Counter()
unique_col_vals = Counter()
protein_keys = Counter()
row_kinds = Counter()
dominant_flag_vals = Counter()

# dominant isoform tracking per protein-ish key
dominant_by_uniprot = defaultdict(list)   # uniprot_id -> list of row identifiers marked dominant
dominant_by_protkey = defaultdict(list)

empty_uniprot = 0
empty_unique = 0
empty_protein_key = 0
empty_ensg = 0
empty_ensp = 0

with open(path, newline='', encoding='utf-8', errors='replace') as f:
    reader = csv.DictReader(f)
    fieldnames = reader.fieldnames
    for row in reader:
        total_rows += 1
        uid = (row.get('uniprot_id') or '').strip()
        uniq = (row.get('UNIQUE') or '').strip()
        pkey = (row.get('protein_key') or '').strip()
        rk = (row.get('row_kind') or '').strip()
        dom = (row.get('dominant_isoform') or '').strip()
        ensg = (row.get('ENSG') or '').strip()
        ensp = (row.get('ENSP') or '').strip()

        if uid:
            uniprot_ids[uid] += 1
        else:
            empty_uniprot += 1

        if uniq:
            unique_col_vals[uniq] += 1
        else:
            empty_unique += 1

        if pkey:
            protein_keys[pkey] += 1
        else:
            empty_protein_key += 1

        if not ensg:
            empty_ensg += 1
        if not ensp:
            empty_ensp += 1

        row_kinds[rk] += 1
        dominant_flag_vals[dom] += 1

        if dom.lower() in ('true', '1', 'yes', 't'):
            if uid:
                dominant_by_uniprot[uid].append(uniq or ensp or 'row')
            if pkey:
                dominant_by_protkey[pkey].append(uniq or ensp or 'row')

print(f"Total data rows: {total_rows}")
print(f"Total columns: {len(fieldnames)}")
print()

print("=== uniprot_id ===")
print(f"Non-empty uniprot_id rows: {total_rows - empty_uniprot}")
print(f"Empty uniprot_id rows: {empty_uniprot}")
print(f"Distinct uniprot_id values: {len(uniprot_ids)}")
dup_uid = {k: v for k, v in uniprot_ids.items() if v > 1}
print(f"uniprot_id values appearing >1 time: {len(dup_uid)} (expected if isoforms share a parent id)")
print()

print("=== UNIQUE column ===")
print(f"Empty UNIQUE rows: {empty_unique}")
print(f"Distinct UNIQUE values: {len(unique_col_vals)}")
dup_unique = {k: v for k, v in unique_col_vals.items() if v > 1}
print(f"UNIQUE values appearing >1 time (should be 0 if truly unique): {len(dup_unique)}")
if dup_unique:
    sample = list(dup_unique.items())[:10]
    print("  sample duplicates:", sample)
print()

print("=== protein_key ===")
print(f"Empty protein_key rows: {empty_protein_key}")
print(f"Distinct protein_key values: {len(protein_keys)}")
print()

print("=== row_kind distribution ===")
for k, v in row_kinds.most_common():
    print(f"  {k!r}: {v}")
print()

print("=== dominant_isoform value distribution ===")
for k, v in dominant_flag_vals.most_common(20):
    print(f"  {k!r}: {v}")
print()

multi_dom_uid = {k: v for k, v in dominant_by_uniprot.items() if len(v) > 1}
print(f"=== uniprot_id groups with >1 row flagged dominant_isoform=True ===")
print(f"Count of such uniprot_ids: {len(multi_dom_uid)}")
if multi_dom_uid:
    for k, v in list(multi_dom_uid.items())[:10]:
        print(f"  {k}: {len(v)} dominant rows -> {v}")
print()

multi_dom_pkey = {k: v for k, v in dominant_by_protkey.items() if len(v) > 1}
print(f"=== protein_key groups with >1 row flagged dominant_isoform=True ===")
print(f"Count of such protein_keys: {len(multi_dom_pkey)}")
if multi_dom_pkey:
    for k, v in list(multi_dom_pkey.items())[:10]:
        print(f"  {k}: {len(v)} dominant rows -> {v}")
print()

total_dominant_true = sum(v for k, v in dominant_flag_vals.items() if k.lower() in ('true','1','yes','t'))
print(f"Total rows with dominant_isoform=True: {total_dominant_true}")
print(f"Total distinct uniprot_id: {len(uniprot_ids)}")
print(f"Total distinct protein_key: {len(protein_keys)}")

print()
print("=== ENSG / ENSP emptiness ===")
print(f"Empty ENSG: {empty_ensg}")
print(f"Empty ENSP: {empty_ensp}")

print()
print("=== Field names ===")
for fn in fieldnames:
    print(" ", fn)
