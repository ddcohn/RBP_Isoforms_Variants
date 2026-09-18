import csv
import sys
import json
import argparse

csv.field_size_limit(sys.maxsize)

TARGET = "/u/project/kappel/ddcohn/RBP/table_260823_with_rna.csv"
GO_RESULT = "/u/home/d/ddcohn/_claude_go_result.json"
TISSUE_RESULT = "/u/home/d/ddcohn/_claude_tissue_result.json"
OT_RESULT = "/u/home/d/ddcohn/_claude_opentargets_result.json"

NEW_COLUMNS = [
    "GO_Cellular_Component",
    "GO_Biological_Process",
    "GO_Molecular_Function",
    "Tissue_Expression_Bulk",
    "OpenTargets_DiseaseId",
    "OpenTargets_DatatypeId",
    "OpenTargets_Score",
]


def fmt_go(pairs):
    return ";".join(f"{gid}:{name}" for gid, name in pairs)


def fmt_tissue(entries):
    parts = []
    for e in entries:
        med = e.get("median")
        med_str = f"{med:.4f}" if isinstance(med, (int, float)) else ""
        parts.append(f"{e.get('tissue','')}|{e.get('datasourceId','')}|{med_str}|{e.get('unit','')}")
    return ";".join(parts)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=TARGET)
    args = ap.parse_args()

    print("Loading GO result...")
    with open(GO_RESULT) as f:
        go_data = json.load(f)
    print(f"  {len(go_data)} proteins")

    print("Loading tissue result...")
    with open(TISSUE_RESULT) as f:
        tissue_data = json.load(f)
    print(f"  {len(tissue_data)} proteins")

    print("Loading OpenTargets result...")
    with open(OT_RESULT) as f:
        ot_data = json.load(f)
    print(f"  {len(ot_data)} proteins")

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

    n_go = n_tissue = n_ot = 0
    for row in rows:
        uid = row["uniprot_accession"]

        go = go_data.get(uid)
        if go:
            row["GO_Cellular_Component"] = fmt_go(go.get("C", []))
            row["GO_Biological_Process"] = fmt_go(go.get("P", []))
            row["GO_Molecular_Function"] = fmt_go(go.get("F", []))
            if go.get("C") or go.get("P") or go.get("F"):
                n_go += 1
        else:
            row["GO_Cellular_Component"] = ""
            row["GO_Biological_Process"] = ""
            row["GO_Molecular_Function"] = ""

        tissue = tissue_data.get(uid)
        if tissue:
            row["Tissue_Expression_Bulk"] = fmt_tissue(tissue)
            if tissue:
                n_tissue += 1
        else:
            row["Tissue_Expression_Bulk"] = ""

        ot = ot_data.get(uid)
        if ot and ot.get("diseaseId"):
            row["OpenTargets_DiseaseId"] = ";".join(ot["diseaseId"])
            row["OpenTargets_DatatypeId"] = ";".join(ot["datatypeId"])
            row["OpenTargets_Score"] = ";".join(ot["score"])
            n_ot += 1
        else:
            row["OpenTargets_DiseaseId"] = ""
            row["OpenTargets_DatatypeId"] = ""
            row["OpenTargets_Score"] = ""

    print(f"\nRows with >=1 GO term: {n_go}/{len(rows)}")
    print(f"Rows with >=1 tissue expression entry: {n_tissue}/{len(rows)}")
    print(f"Rows with >=1 OpenTargets disease association: {n_ot}/{len(rows)}")

    print(f"\nWriting merged table to {args.out} ...")
    with open(args.out, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=out_fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    print("Done.")


if __name__ == "__main__":
    main()
