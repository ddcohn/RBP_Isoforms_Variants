import csv
import sys
import gzip

csv.field_size_limit(sys.maxsize)

SRC = "/u/project/kappel/ddcohn/variant_summary.txt.gz"
OUT = "/u/project/kappel/ddcohn/ClinVar_variant_summary_GRCh38.csv"

KEEP_COLS = [
    "VariationID", "AlleleID", "Type", "Name",
    "GeneID", "GeneSymbol", "HGNC_ID",
    "ClinicalSignificance", "ClinSigSimple", "ReviewStatus", "NumberSubmitters",
    "PhenotypeList", "Origin", "OriginSimple",
    "Chromosome", "PositionVCF", "ReferenceAlleleVCF", "AlternateAlleleVCF",
    "SomaticClinicalImpact", "Oncogenicity",
]

n_total = 0
n_grch38 = 0
seen_variation_ids = set()
n_dup_variation_id = 0

with gzip.open(SRC, mode="rt", newline="", encoding="utf-8", errors="replace") as f, \
     open(OUT, "w", newline="", encoding="utf-8") as out:
    reader = csv.DictReader(f, delimiter="\t")
    # the header starts with '#AlleleID' - fix the key
    reader.fieldnames = [fn.lstrip("#") for fn in reader.fieldnames]
    writer = csv.DictWriter(out, fieldnames=KEEP_COLS)
    writer.writeheader()
    for row in reader:
        n_total += 1
        if row.get("Assembly") != "GRCh38":
            continue
        n_grch38 += 1
        vid = row.get("VariationID", "")
        if vid in seen_variation_ids:
            n_dup_variation_id += 1
        else:
            seen_variation_ids.add(vid)
        writer.writerow({c: row.get(c, "") for c in KEEP_COLS})
        if n_total % 2_000_000 == 0:
            print(f"  ... {n_total:,} total rows scanned, {n_grch38:,} GRCh38 rows written", flush=True)

print(f"\nTotal rows scanned: {n_total:,}")
print(f"GRCh38 rows written: {n_grch38:,}")
print(f"Distinct VariationID among GRCh38 rows: {len(seen_variation_ids):,}")
print(f"Duplicate VariationID rows within GRCh38 subset: {n_dup_variation_id:,}")
