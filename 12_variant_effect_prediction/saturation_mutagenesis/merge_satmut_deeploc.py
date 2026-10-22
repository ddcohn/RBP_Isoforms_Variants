"""
Concatenates all 226 satmut_deeploc_results/chunk_NNNN/results_*.csv
files (each already verified complete via a full-array content audit)
into one single mutation-level DeepLoc results file. Re-parsed and
re-written as TSV (source is CSV; some fields contain embedded "|"
lists but values were checked to never contain literal tabs, so TSV
output is safe and consistent with the other three satmut result
files for the later 4-way join by mutation_id).
"""
import csv
import glob
import sys

csv.field_size_limit(sys.maxsize)

RESULTS_DIR = "/u/project/kappel/ddcohn/protein_variant_effects/rbp_only/satmut_deeploc_results"
OUT_FILE = "/u/project/kappel/ddcohn/protein_variant_effects/rbp_only/satmut_deeploc_all.tsv"

chunk_dirs = sorted(glob.glob(f"{RESULTS_DIR}/chunk_*"))
print(f"Found {len(chunk_dirs)} chunk dirs")
assert len(chunk_dirs) == 226, f"expected 226 chunks, found {len(chunk_dirs)}"

n_written = 0
header_written = False
with open(OUT_FILE, "w", newline="") as out:
    writer = csv.writer(out, delimiter="\t")
    for d in chunk_dirs:
        matches = sorted(glob.glob(f"{d}/results_*.csv"))
        assert len(matches) == 1, f"expected exactly 1 results file in {d}, found {len(matches)}"
        fp = matches[0]
        with open(fp, newline="") as f:
            reader = csv.reader(f)
            header = next(reader)
            if not header_written:
                writer.writerow(header)
                header_written = True
            for row in reader:
                assert not any("\t" in v for v in row), f"literal tab found in a field in {fp}"
                writer.writerow(row)
                n_written += 1

print(f"Rows written: {n_written}")
print(f"Wrote {OUT_FILE}")
