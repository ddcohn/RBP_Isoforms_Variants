import csv
import sys
import os

csv.field_size_limit(sys.maxsize)

# usage: build_clinvar_spliceai_chunks.py <raw_csv> <ref_fasta> <outdir> <n_chunks>
SRC = sys.argv[1]
FASTA = sys.argv[2]
OUTDIR = sys.argv[3]
N_CHUNKS = int(sys.argv[4])

os.makedirs(OUTDIR, exist_ok=True)

CHROM_MAP = {"MT": "M"}
VALID_CHROMS = set([str(i) for i in range(1, 23)] + ["X", "Y", "MT"])

# get lengths for all primary contigs
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

rows = []
skipped_un = 0
skipped_bad_pos = 0
with open(SRC, newline="", encoding="utf-8") as f:
    for row in csv.DictReader(f):
        if row["CoordinateAssembly"] != "GRCh38":
            continue
        chrom_raw = row["Chromosome"]
        if chrom_raw not in VALID_CHROMS:
            skipped_un += 1
            continue
        chrom = "chr" + CHROM_MAP.get(chrom_raw, chrom_raw)
        pos = row["PositionVCF"]
        try:
            if int(pos) < 1:
                skipped_bad_pos += 1
                continue
        except ValueError:
            skipped_bad_pos += 1
            continue
        vid = row["VariationID"]
        ref = row["ReferenceAlleleVCF"]
        alt = row["AlternateAlleleVCF"]
        rows.append((chrom, pos, vid, ref, alt))

print(f"Total usable GRCh38 variants: {len(rows)} (skipped {skipped_un} unplaced, "
      f"skipped {skipped_bad_pos} with invalid PositionVCF)")

n = len(rows)
chunk_size = (n + N_CHUNKS - 1) // N_CHUNKS

for i in range(N_CHUNKS):
    start = i * chunk_size
    end = min(start + chunk_size, n)
    if start >= n:
        break
    chunk_rows = rows[start:end]
    out_path = os.path.join(OUTDIR, f"chunk_{i+1:03d}.vcf")
    with open(out_path, "w") as out:
        out.writelines(header_lines)
        for chrom, pos, vid, ref, alt in chunk_rows:
            out.write(f"{chrom}\t{pos}\t{vid}\t{ref}\t{alt}\t.\t.\t.\n")
    print(f"Wrote {out_path}: {len(chunk_rows)} variants")

print("Done.")
