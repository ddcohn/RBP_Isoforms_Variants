import csv
import glob
import sys

csv.field_size_limit(sys.maxsize)

RESULTS_DIR = "/u/project/kappel/ddcohn/protein_variant_effects/cmc_protgps_results"
MAP_FILE = "/u/project/kappel/ddcohn/protein_variant_effects/cmc_variant_sequence_map.tsv"
OUT_FILE = "/u/project/kappel/ddcohn/protein_variant_effects/cmc_protgps_deltas.tsv"

COMPARTMENTS = [
    "NUCLEAR_SPECKLE", "P-BODY", "PML-BDOY", "POST_SYNAPTIC_DENSITY",
    "STRESS_GRANULE", "CHROMOSOME", "NUCLEOLUS", "NUCLEAR_PORE_COMPLEX",
    "CAJAL_BODY", "RNA_GRANULE", "CELL_JUNCTION", "TRANSCRIPTIONAL",
]

print("Loading all CMC protGPS result chunks...")
scores = {}
files = sorted(glob.glob(f"{RESULTS_DIR}/chunk_*.tsv"))
for fp in files:
    with open(fp) as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader)
        for row in reader:
            scores[row[0]] = row[1:]

print(f"Loaded {len(scores)} sequence predictions")

n_written = 0
n_missing_wt = 0
n_missing_mut = 0
n_na = 0

with open(MAP_FILE, newline="") as f, open(OUT_FILE, "w", newline="") as out:
    reader = csv.DictReader(f, delimiter="\t")
    writer = csv.writer(out, delimiter="\t")
    header = ["GENOMIC_MUTATION_ID", "GENE_NAME", "accession", "category"]
    header += [f"WT_{c}" for c in COMPARTMENTS]
    header += [f"MUT_{c}" for c in COMPARTMENTS]
    header += [f"DELTA_{c}" for c in COMPARTMENTS]
    writer.writerow(header)

    for row in reader:
        wt_scores = scores.get(row["accession"])
        mut_scores = scores.get(row["mutant_id"])
        if wt_scores is None:
            n_missing_wt += 1
            continue
        if mut_scores is None:
            n_missing_mut += 1
            continue
        if "NA" in wt_scores or "NA" in mut_scores:
            n_na += 1
            continue
        wt_f = [float(v) for v in wt_scores]
        mut_f = [float(v) for v in mut_scores]
        delta = [m - w for w, m in zip(wt_f, mut_f)]
        out_row = [row["GENOMIC_MUTATION_ID"], row["GENE_NAME"], row["accession"], row["category"]]
        out_row += [f"{v:.4f}" for v in wt_f]
        out_row += [f"{v:.4f}" for v in mut_f]
        out_row += [f"{v:.4f}" for v in delta]
        writer.writerow(out_row)
        n_written += 1

print(f"\nRows written: {n_written}")
print(f"Missing WT score: {n_missing_wt}")
print(f"Missing mutant score: {n_missing_mut}")
print(f"Skipped (NA, too-long sequence): {n_na}")
print(f"\nWrote {OUT_FILE}")
