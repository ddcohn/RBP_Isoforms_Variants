import csv
import sys

csv.field_size_limit(sys.maxsize)

# usage: extract_clinvar_metadata.py <raw_csv> <out_tsv>
SRC = sys.argv[1]
OUT = sys.argv[2]

n = 0
with open(SRC, newline="", encoding="utf-8") as f, open(OUT, "w", newline="") as out:
    reader = csv.DictReader(f)
    writer = csv.writer(out, delimiter="\t")
    writer.writerow(["VariationID", "GeneSymbol", "ClinicalSignificance"])
    for row in reader:
        if row["CoordinateAssembly"] != "GRCh38":
            continue
        writer.writerow([row["VariationID"], row["GeneSymbol"], row["ClinicalSignificance"]])
        n += 1

print(f"Wrote {n} rows to {OUT}")
