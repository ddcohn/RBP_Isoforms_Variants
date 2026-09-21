import csv
import sys
import gzip
from collections import Counter

csv.field_size_limit(sys.maxsize)

SRC = "/u/project/kappel/ddcohn/Cosmic_GenomeScreensMutant_v104_GRCh38.tsv.gz"

seen = set()
n_rows = 0
c = Counter()

with gzip.open(SRC, "rt", encoding="utf-8", errors="replace") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        n_rows += 1
        chrom = row.get("CHROMOSOME", "")
        pos = row.get("GENOME_START", "")
        ref = row.get("GENOMIC_WT_ALLELE", "")
        alt = row.get("GENOMIC_MUT_ALLELE", "")
        if not chrom or not pos or not ref or not alt or ref == "-" or alt == "-":
            c["unusable"] += 1
            continue
        key = (chrom, pos, ref, alt)
        seen.add(key)
        if n_rows % 5000000 == 0:
            print(f"  processed {n_rows} rows, {len(seen)} unique so far...")

print(f"\nTotal rows: {n_rows}")
print(f"Unusable (missing/'-' alleles): {c['unusable']}")
print(f"Unique (chrom,pos,ref,alt) variants: {len(seen)}")
