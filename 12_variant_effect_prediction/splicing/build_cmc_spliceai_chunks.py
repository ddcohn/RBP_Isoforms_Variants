import csv
import sys
import gzip
import os
import re

csv.field_size_limit(sys.maxsize)

SRC = "/u/project/kappel/ddcohn/CancerMutationCensus_AllData_v104_GRCh37.tsv.gz"
FASTA = "/u/project/kappel/ddcohn/RNA-GPS/rnagps/reference/GRCh38.primary_assembly.genome.fa"
OUTDIR = "/u/project/kappel/ddcohn/SpliceAI/cmc_run/chunks"
N_CHUNKS = 260

os.makedirs(OUTDIR, exist_ok=True)

VALID_CHROMS = set([str(i) for i in range(1, 23)] + ["X", "Y", "MT"])
CHROM_MAP = {"MT": "M"}
POS_RE = re.compile(r"^([0-9XYMT]+):(\d+)-(\d+)$")

# contig lengths
lengths = {}
name = None
length = 0
with open(FASTA) as f:
    for line in f:
        if line.startswith(">"):
            if name is not None:
                lengths[name] = length
            name = line[1:].split()[0]
            length = 0
        else:
            length += len(line.strip())
    if name is not None:
        lengths[name] = length

contig_order = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"]
header_lines = ["##fileformat=VCFv4.2\n"]
for c in contig_order:
    header_lines.append(f"##contig=<ID={c},length={lengths[c]}>\n")
header_lines.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

seen = set()
rows = []
n_total = 0
skipped_no_pos = 0
skipped_bad_allele = 0
skipped_dup = 0

with gzip.open(SRC, "rt", encoding="utf-8", errors="replace") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        n_total += 1
        pos38 = row.get("Mutation genome position GRCh38", "")
        if not pos38:
            skipped_no_pos += 1
            continue
        m = POS_RE.match(pos38)
        if not m:
            skipped_no_pos += 1
            continue
        chrom_raw, start, stop = m.group(1), m.group(2), m.group(3)
        if chrom_raw not in VALID_CHROMS:
            skipped_no_pos += 1
            continue
        ref = row.get("GENOMIC_WT_ALLELE_SEQ", "")
        alt = row.get("GENOMIC_MUT_ALLELE_SEQ", "")
        if not ref or not alt or ref == "-" or alt == "-":
            skipped_bad_allele += 1
            continue
        chrom = "chr" + CHROM_MAP.get(chrom_raw, chrom_raw)
        key = (chrom, start, ref, alt)
        if key in seen:
            skipped_dup += 1
            continue
        seen.add(key)
        mid = row.get("GENOMIC_MUTATION_ID", "NA")
        rows.append((chrom, start, mid, ref, alt))
        if n_total % 1000000 == 0:
            print(f"  processed {n_total} rows, {len(rows)} unique usable so far...")

print(f"\nTotal rows: {n_total}")
print(f"Skipped (no/bad GRCh38 pos): {skipped_no_pos}")
print(f"Skipped (bad allele): {skipped_bad_allele}")
print(f"Skipped (duplicate variant): {skipped_dup}")
print(f"Usable unique variants: {len(rows)}")

n = len(rows)
chunk_size = (n + N_CHUNKS - 1) // N_CHUNKS

for i in range(N_CHUNKS):
    start_i = i * chunk_size
    end_i = min(start_i + chunk_size, n)
    if start_i >= n:
        break
    chunk_rows = rows[start_i:end_i]
    out_path = os.path.join(OUTDIR, f"chunk_{i+1:03d}.vcf")
    with open(out_path, "w") as out:
        out.writelines(header_lines)
        for chrom, pos, mid, ref, alt in chunk_rows:
            out.write(f"{chrom}\t{pos}\t{mid}\t{ref}\t{alt}\t.\t.\t.\n")
    print(f"Wrote {out_path}: {len(chunk_rows)} variants")

print("Done.")
