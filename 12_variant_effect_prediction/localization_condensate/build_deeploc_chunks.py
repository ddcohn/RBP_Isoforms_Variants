import os
import sys

# usage: build_deeploc_chunks.py <wt_fasta> <mutant_fasta> <outdir> <n_chunks>
# shared by ClinVar and CMC -- takes whichever WT/mutant FASTA pair is passed in.
WT_FASTA = sys.argv[1]
MUT_FASTA = sys.argv[2]
OUTDIR = sys.argv[3]
N_CHUNKS = int(sys.argv[4])

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
