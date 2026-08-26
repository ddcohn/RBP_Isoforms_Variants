import csv
import sys
import json

csv.field_size_limit(sys.maxsize)

TABLE = "/u/project/kappel/ddcohn/RBP/table_260823_with_rna.csv"

with open("/u/home/d/ddcohn/ppi_local_v125_resolved_result.json") as f:
    local_result = json.load(f)  # uid -> [partner_uids, in_df_uids]

with open("/u/home/d/ddcohn/ppi_api_v125_result.json") as f:
    api_result = json.load(f)

rows = []
with open(TABLE, newline="", encoding="utf-8", errors="replace") as f:
    reader = csv.DictReader(f)
    fieldnames = reader.fieldnames
    for row in reader:
        rows.append(row)

new_cols = [
    "PPI_UniProt_Partners_v125_local",
    "PPI_UniProt_Partners_in_Dataframe_v125_local",
    "PPI_UniProt_Partners_v125_api",
    "PPI_UniProt_Partners_in_Dataframe_v125_api",
]
out_fields = fieldnames + [c for c in new_cols if c not in fieldnames]

for row in rows:
    uid = row["uniprot_accession"].strip()
    l_partners, l_df = local_result.get(uid, ([], []))
    a_partners, a_df = api_result.get(uid, ([], []))
    row["PPI_UniProt_Partners_v125_local"] = ";".join(l_partners)
    row["PPI_UniProt_Partners_in_Dataframe_v125_local"] = ";".join(l_df)
    row["PPI_UniProt_Partners_v125_api"] = ";".join(a_partners)
    row["PPI_UniProt_Partners_in_Dataframe_v125_api"] = ";".join(a_df)

with open(TABLE, "w", newline="", encoding="utf-8") as f:
    writer = csv.DictWriter(f, fieldnames=out_fields)
    writer.writeheader()
    writer.writerows(rows)

n_local = sum(1 for r in rows if r["PPI_UniProt_Partners_v125_local"])
n_api = sum(1 for r in rows if r["PPI_UniProt_Partners_v125_api"])
n_v120 = sum(1 for r in rows if r.get("PPI_UniProt_Partners"))
print(f"Total rows: {len(rows)}")
print(f"v12.0 (existing, live API): {n_v120} with partners")
print(f"v12.5 local file:            {n_local} with partners")
print(f"v12.5 live API:               {n_api} with partners")
