import csv
import json

RESULT = "/u/home/d/ddcohn/_claude_go_result.json"
OUT = "/u/home/d/ddcohn/_claude_go_table.csv"

with open(RESULT) as f:
    data = json.load(f)


def fmt(pairs):
    return ";".join(f"{gid}:{name}" for gid, name in pairs)


with open(OUT, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["uniprot_accession", "GO_Cellular_Component", "GO_Biological_Process", "GO_Molecular_Function"])
    for acc, v in data.items():
        writer.writerow([acc, fmt(v.get("C", [])), fmt(v.get("P", [])), fmt(v.get("F", []))])

print(f"Wrote {len(data)} rows to {OUT}")
