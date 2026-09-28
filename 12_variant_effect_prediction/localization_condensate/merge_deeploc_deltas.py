import csv
import glob
import os
import sys

# usage: merge_deeploc_deltas.py clinvar|cmc <base_dir>
which = sys.argv[1]
BASE = sys.argv[2]
if which == "clinvar":
    RESULTS, CHUNKS = f"{BASE}/deeploc_results", f"{BASE}/deeploc_chunks"
    MAP_FILE, OUT_FILE = f"{BASE}/variant_sequence_map.tsv", f"{BASE}/clinvar_deeploc_deltas.tsv"
    ID_COLS = ["VariationID", "GeneID", "GeneSymbol", "accession", "category"]
elif which == "cmc":
    RESULTS, CHUNKS = f"{BASE}/cmc_deeploc_results", f"{BASE}/cmc_deeploc_chunks"
    MAP_FILE, OUT_FILE = f"{BASE}/cmc_variant_sequence_map.tsv", f"{BASE}/cmc_deeploc_deltas.tsv"
    ID_COLS = ["GENOMIC_MUTATION_ID", "GENE_NAME", "accession", "category"]
else:
    sys.exit("argument must be clinvar or cmc")

csv.field_size_limit(sys.maxsize)
LOCS = ["Cytoplasm", "Nucleus", "Extracellular", "Cell membrane", "Mitochondrion", "Plastid",
        "Endoplasmic reticulum", "Lysosome/Vacuole", "Golgi apparatus", "Peroxisome"]
MEMB = ["Peripheral", "Transmembrane", "Lipid anchor", "Soluble"]
NUM = LOCS + MEMB


def fasta_ids(path):
    with open(path) as f:
        return [l[1:].strip() for l in f if l.startswith(">")]


print("Loading DeepLoc results and verifying each chunk against its input FASTA...")
scores, labels = {}, {}
n_chunks = n_problem = 0
for cdir in sorted(glob.glob(f"{RESULTS}/chunk_*")):
    name = os.path.basename(cdir)
    csvs = glob.glob(f"{cdir}/results_*.csv")
    ids = fasta_ids(f"{CHUNKS}/{name}.fasta")
    got = set()
    for fp in csvs:
        with open(fp, newline="") as f:
            for row in csv.DictReader(f):
                pid = row["Protein_ID"]
                got.add(pid)
                scores[pid] = [float(row[c]) for c in NUM]
                labels[pid] = row["Localizations"]
    missing = set(ids) - got
    dup = len(ids) - len(set(ids))
    n_chunks += 1
    if len(csvs) != 1 or missing:
        n_problem += 1
        print(f"  PROBLEM {name}: csv_files={len(csvs)} missing_ids={len(missing)}")
    elif dup:
        print(f"  note {name}: {dup} duplicated FASTA ID(s); DeepLoc keeps one row per ID")
print(f"Chunks checked: {n_chunks}, with problems: {n_problem}")
if n_problem:
    sys.exit("Aborting: fix incomplete chunks first")
print(f"Loaded {len(scores)} sequence predictions")

n_written = n_missing_wt = n_missing_mut = 0
with open(MAP_FILE, newline="") as f, open(OUT_FILE, "w", newline="") as out:
    reader = csv.DictReader(f, delimiter="\t")
    w = csv.writer(out, delimiter="\t")
    header = ID_COLS + ["WT_label", "MUT_label", "label_changed"]
    header += [f"WT_{c}" for c in NUM] + [f"MUT_{c}" for c in NUM] + [f"DELTA_{c}" for c in NUM]
    w.writerow(header)
    for row in reader:
        wt, mut = scores.get(row["accession"]), scores.get(row["mutant_id"])
        if wt is None:
            n_missing_wt += 1
            continue
        if mut is None:
            n_missing_mut += 1
            continue
        wl, ml = labels[row["accession"]], labels[row["mutant_id"]]
        out_row = [row[c] for c in ID_COLS] + [wl, ml, int(wl != ml)]
        out_row += [f"{v:.4f}" for v in wt] + [f"{v:.4f}" for v in mut]
        out_row += [f"{m - x:.4f}" for x, m in zip(wt, mut)]
        w.writerow(out_row)
        n_written += 1

print(f"Rows written: {n_written}")
print(f"Missing WT prediction: {n_missing_wt}")
print(f"Missing mutant prediction: {n_missing_mut}")
print(f"Wrote {OUT_FILE}")
