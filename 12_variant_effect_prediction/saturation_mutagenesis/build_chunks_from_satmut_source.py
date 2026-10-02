import csv
import gzip
import os
import sys

csv.field_size_limit(sys.maxsize)

SRC = "/u/project/kappel/asharma/RBP/saturation_mutagenesis/rbp_satmut_all.csv.gz"
OUTDIR = "/u/project/kappel/ddcohn/protein_variant_effects/rbp_only/satmut_chunks"
MANIFEST = "/u/project/kappel/ddcohn/protein_variant_effects/rbp_only/satmut_manifest.tsv"
os.makedirs(OUTDIR, exist_ok=True)

TARGET_CHUNK_SIZE = 140_000


def open_chunk(i):
    return open(os.path.join(OUTDIR, f"chunk_{i:04d}.fasta"), "w")


chunk_idx = 1
chunk_n = 0
out = open_chunk(chunk_idx)
n_total = 0
n_written = 0
n_wt_skipped = 0
n_proteins = set()

manifest = open(MANIFEST, "w")
manifest.write("uniprot_accession\tRBP_name\tn_mutants\n")
protein_counts = {}

with gzip.open(SRC, "rt", encoding="utf-8", errors="replace") as f:
    r = csv.DictReader(f)
    for row in r:
        n_total += 1
        if row["is_wt"] == "True":
            n_wt_skipped += 1
            continue
        acc = row["uniprot_accession"]
        gene = row["RBP_name"]
        mut = row["mutation"]
        seq = row["sequence"]
        n_proteins.add(acc)
        protein_counts[(acc, gene)] = protein_counts.get((acc, gene), 0) + 1
        mut_id = f"{acc}|{gene}|{mut}"
        out.write(f">{mut_id}\n{seq}\n")
        chunk_n += 1
        n_written += 1
        if chunk_n >= TARGET_CHUNK_SIZE:
            out.close()
            chunk_idx += 1
            chunk_n = 0
            out = open_chunk(chunk_idx)
        if n_total % 2_000_000 == 0:
            print(f"  processed {n_total:,} rows, written {n_written:,} mutant sequences so far...")

out.close()

for (acc, gene), c in sorted(protein_counts.items()):
    manifest.write(f"{acc}\t{gene}\t{c}\n")
manifest.close()

print(f"\nTotal rows in source file: {n_total:,}")
print(f"WT reference rows skipped: {n_wt_skipped:,}")
print(f"Mutant sequences written: {n_written:,}")
print(f"Unique proteins: {len(n_proteins)}")
print(f"Chunks: {chunk_idx}")
print(f"Wrote {MANIFEST}")
