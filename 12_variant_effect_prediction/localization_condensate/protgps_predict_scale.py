import sys
import os
import pickle
import time
from argparse import Namespace

import torch

PROTGPS_DIR = "/u/project/kappel/ddcohn/localization_tools/protgps"
sys.path.append(PROTGPS_DIR)
from protgps.utils.loading import get_object

# usage: protgps_predict_scale.py <chunk_id> <in_fasta> <out_tsv>
# shared by ClinVar and CMC -- takes whichever chunk FASTA is passed in.
CHUNK = sys.argv[1]
IN_FASTA = sys.argv[2]
OUT_TSV = sys.argv[3]
os.makedirs(os.path.dirname(OUT_TSV), exist_ok=True)

device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
print("Using device:", device)

COMPARTMENT_CLASSES = [
    "nuclear_speckle", "p-body", "pml-bdoy", "post_synaptic_density",
    "stress_granule", "chromosome", "nucleolus", "nuclear_pore_complex",
    "cajal_body", "rna_granule", "cell_junction", "transcriptional",
]

MAX_LEN = 5000  # ESM2 8M / attention memory guard; longer sequences skipped
BASE_BATCH = 8


def load_model(snargs):
    model = get_object(snargs.lightning_name, "lightning")(snargs)
    model = model.load_from_checkpoint(
        checkpoint_path=snargs.model_path,
        strict=not snargs.relax_checkpoint_matching,
        **{"args": snargs},
    )
    return model


@torch.no_grad()
def predict_batch(model, sequences):
    out = model.model({"x": sequences})
    return torch.sigmoid(out["logit"]).to("cpu")


args = Namespace(**pickle.load(open(
    os.path.join(PROTGPS_DIR, "checkpoints/protgps/32bf44b16a4e770a674896b81dfb3729.args"), "rb")))
args.model_path = os.path.join(
    PROTGPS_DIR, "checkpoints/protgps/32bf44b16a4e770a674896b81dfb3729epoch=26.ckpt")
args.pretrained_hub_dir = "/u/project/kappel/ddcohn/model_cache/torch_hub"

model = load_model(args)
model.eval()
model = model.to(device)

# load sequences
records = []
with open(IN_FASTA) as f:
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

print(f"Loaded {len(records)} sequences")

skipped = [r for r in records if len(r[1]) > MAX_LEN]
usable = [r for r in records if len(r[1]) <= MAX_LEN]
print(f"Skipped (>{MAX_LEN} aa): {len(skipped)}; usable: {len(usable)}")

# sort by length so batches have similar padding needs
usable.sort(key=lambda r: len(r[1]))

start = time.time()
with open(OUT_TSV, "w") as out:
    out.write("id\t" + "\t".join(c.upper() + "_Score" for c in COMPARTMENT_CLASSES) + "\n")
    for skip_id, skip_seq in skipped:
        out.write(skip_id + "\t" + "\t".join(["NA"] * len(COMPARTMENT_CLASSES)) + "\n")

    i = 0
    n = len(usable)
    while i < n:
        # adapt batch size down for longer sequences to control memory
        cur_len = len(usable[i][1])
        if cur_len > 2000:
            bs = 2
        elif cur_len > 1000:
            bs = 4
        else:
            bs = BASE_BATCH
        batch = usable[i:i + bs]
        seqs = [s for _, s in batch]
        try:
            scores = predict_batch(model, seqs)
        except torch.cuda.OutOfMemoryError:
            torch.cuda.empty_cache()
            scores = None
            for single_id, single_seq in batch:
                try:
                    s = predict_batch(model, [single_seq])
                    out.write(single_id + "\t" + "\t".join(f"{v:.4f}" for v in s[0].tolist()) + "\n")
                except torch.cuda.OutOfMemoryError:
                    torch.cuda.empty_cache()
                    out.write(single_id + "\t" + "\t".join(["OOM"] * len(COMPARTMENT_CLASSES)) + "\n")
            i += bs
            continue
        for (rid, _), row in zip(batch, scores):
            out.write(rid + "\t" + "\t".join(f"{v:.4f}" for v in row.tolist()) + "\n")
        i += bs
        if i % 5000 < bs:
            elapsed = time.time() - start
            print(f"  processed {i}/{n} usable, elapsed {elapsed:.0f}s, rate {i/elapsed:.1f} seqs/sec")

elapsed = time.time() - start
print(f"Done. {n} usable sequences in {elapsed:.0f}s ({n/elapsed:.1f} seqs/sec)")
print(f"Wrote {OUT_TSV}")
