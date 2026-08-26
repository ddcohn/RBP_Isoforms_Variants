import csv
import sys
from collections import defaultdict

path = "/u/project/kappel/RBP/Isoform_Table/Isoform_Post_Merge_PSLab_OpenTargets_Updated_Interim_20260817.csv"
csv.field_size_limit(sys.maxsize)

rows_per_uid = defaultdict(int)
dominant_per_uid = defaultdict(int)
gene_per_uid = {}

with open(path, newline='', encoding='utf-8', errors='replace') as f:
    reader = csv.DictReader(f)
    for row in reader:
        uid = (row.get('uniprot_id') or '').strip()
        if not uid:
            continue
        rows_per_uid[uid] += 1
        dom = (row.get('dominant_isoform') or '').strip()
        if dom in ('1', 'True', 'true'):
            dominant_per_uid[uid] += 1
        gs = (row.get('gene_symbol') or '').strip()
        if gs and uid not in gene_per_uid:
            gene_per_uid[uid] = gs

no_dominant = [uid for uid in rows_per_uid if dominant_per_uid.get(uid, 0) == 0]

print(f"Total distinct uniprot_id: {len(rows_per_uid)}")
print(f"uniprot_id with ZERO rows marked dominant_isoform=1: {len(no_dominant)}")
print()
print("Sample (first 25) of proteins with no dominant isoform assigned:")
for uid in no_dominant[:25]:
    print(f"  {uid}  gene={gene_per_uid.get(uid,'?')}  total_rows={rows_per_uid[uid]}")

# distribution of how many isoform rows these "no dominant" proteins have
from collections import Counter
dist = Counter(rows_per_uid[uid] for uid in no_dominant)
print()
print("Row-count distribution among 'no dominant' proteins (n_rows -> count of proteins):")
for k in sorted(dist):
    print(f"  {k} rows: {dist[k]} proteins")
