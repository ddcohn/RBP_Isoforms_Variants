import os

WT_FASTA = "/u/project/kappel/ddcohn/protein_variant_effects/cmc_wt_sequences.fasta"
MUT_FASTA = "/u/project/kappel/ddcohn/protein_variant_effects/cmc_mutant_sequences.fasta"
OUTDIR = "/u/project/kappel/ddcohn/protein_variant_effects/cmc_deeploc_chunks"
N_CHUNKS = 32

os.makedirs(OUTDIR, exist_ok=True)


def read_fasta(path):
    records = []
    with open(path) as f:
        header = None
        seq = []
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq)))
                header = line[1:]
                seq = []
            else:
                seq.append(line)
        if header is not None:
            records.append((header, "".join(seq)))
    return records


records = read_fasta(WT_FASTA) + read_fasta(MUT_FASTA)
print(f"Total sequences: {len(records)}")

n = len(records)
chunk_size = (n + N_CHUNKS - 1) // N_CHUNKS

for i in range(N_CHUNKS):
    start = i * chunk_size
    end = min(start + chunk_size, n)
    if start >= n:
        break
    chunk = records[start:end]
    out_path = os.path.join(OUTDIR, f"chunk_{i+1:03d}.fasta")
    with open(out_path, "w") as out:
        for header, seq in chunk:
            out.write(f">{header}\n{seq}\n")
    print(f"Wrote {out_path}: {len(chunk)} sequences")

print("Done.")
