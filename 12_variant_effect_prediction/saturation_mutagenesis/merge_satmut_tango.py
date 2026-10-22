"""
Concatenates all 226 satmut_tango_results/chunk_NNNN.tsv files (each
already verified complete: real row count == real input FASTA count,
per-chunk) into one single mutation-level TANGO results file.

No WT-vs-MUT delta computed here -- that's a separate step (would need
a WT TANGO baseline per protein, which hasn't been run). This step is
just the straight concatenation the saturation_mutagenesis README's
"Merging results" section describes.
"""
import csv
import glob
import sys

csv.field_size_limit(sys.maxsize)

RESULTS_DIR = "/u/project/kappel/ddcohn/protein_variant_effects/rbp_only/satmut_tango_results"
OUT_FILE = "/u/project/kappel/ddcohn/protein_variant_effects/rbp_only/satmut_tango_all.tsv"

files = sorted(glob.glob(f"{RESULTS_DIR}/chunk_*.tsv"))
print(f"Found {len(files)} chunk files")
assert len(files) == 226, f"expected 226 chunks, found {len(files)}"

n_written = 0
header_written = False
with open(OUT_FILE, "w", newline="") as out:
    writer = csv.writer(out, delimiter="\t")
    for fp in files:
        with open(fp, newline="") as f:
            reader = csv.reader(f, delimiter="\t")
            header = next(reader)
            if not header_written:
                writer.writerow(header)
                header_written = True
            for row in reader:
                writer.writerow(row)
                n_written += 1

print(f"Rows written: {n_written}")
print(f"Wrote {OUT_FILE}")
