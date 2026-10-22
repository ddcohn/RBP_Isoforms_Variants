"""
Concatenates all 226 satmut_protgps_results/chunk_NNNN.tsv files (each
already verified complete via a full-array content audit: real output
row count == real input FASTA count, per chunk) into one single
mutation-level protGPS results file. Straight concatenation, in chunk
order -- protgps_predict_scale.py writes each chunk fresh in one pass
(no resume/append logic like TANGO's), so within-chunk row order
already matches the original FASTA order.
"""
import csv
import glob
import sys

csv.field_size_limit(sys.maxsize)

RESULTS_DIR = "/u/project/kappel/ddcohn/protein_variant_effects/rbp_only/satmut_protgps_results"
OUT_FILE = "/u/project/kappel/ddcohn/protein_variant_effects/rbp_only/satmut_protgps_all.tsv"

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
