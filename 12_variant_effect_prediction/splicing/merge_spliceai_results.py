import csv
import glob

RESULTS_DIR = "/u/project/kappel/ddcohn/SpliceAI/full_run/results"
OUT = "/u/project/kappel/ddcohn/protein_variant_effects/clinvar_spliceai_scores.tsv"

SPLICEAI_COLS = ["DS_AG", "DS_AL", "DS_DG", "DS_DL", "DP_AG", "DP_AL", "DP_DG", "DP_DL"]

n_lines = 0
n_missing = 0

files = sorted(glob.glob(f"{RESULTS_DIR}/chunk_*_out.vcf"))
print(f"Found {len(files)} result files")

with open(OUT, "w", newline="") as out:
    writer = csv.writer(out, delimiter="\t")
    writer.writerow(["VariationID"] + SPLICEAI_COLS)
    for fp in files:
        with open(fp) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                n_lines += 1
                fields = line.rstrip("\n").split("\t")
                vid = fields[2]
                info = fields[7]
                if not info.startswith("SpliceAI="):
                    n_missing += 1
                    continue
                parts = info[len("SpliceAI="):].split("|")
                if len(parts) != 10:
                    n_missing += 1
                    continue
                vals = parts[2:10]
                if all(v == "." for v in vals):
                    n_missing += 1
                    continue
                writer.writerow([vid] + vals)

print(f"Total variant lines processed: {n_lines}")
print(f"Missing/unscored: {n_missing}")
print(f"Scored variants written: {n_lines - n_missing}")
print(f"Wrote {OUT}")
