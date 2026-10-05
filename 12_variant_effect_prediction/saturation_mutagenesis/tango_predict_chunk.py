"""
Run TANGO (aggregation propensity) on every sequence in one satmut chunk
FASTA, writing one summary row per sequence to a TSV.

Resumable by design: TANGO's per-call cost scales worse than linearly
with sequence length (measured ~0.35s at 594 residues, ~7.3s at 8797
residues), and chunks group mutants by protein, so a chunk dominated by
one very long protein's mutants could in the worst case take far longer
than a typical chunk. Rather than size h_rt for that worst case, this
script checks which mutation_ids are already present in the output TSV
on startup and skips them -- so if a task is killed by its wall-time
limit partway through, simply resubmitting the same task continues
where it left off instead of losing the whole chunk (unlike the
DeepLoc failure mode documented in README.md, which required a separate
backfill job because DeepLoc batches internally and has no equivalent
resume point).

Usage: python3 tango_predict_chunk.py <chunk_fasta> <out_tsv> <tango_binary>
"""
import os
import subprocess
import sys
import uuid

TANGO_FIELDS = ["AGG", "AMYLO", "TURN", "HELIX", "HELAGG", "BETA"]


def read_fasta(path):
    header = None
    seq_parts = []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_parts)
                header = line[1:]
                seq_parts = []
            else:
                seq_parts.append(line)
    if header is not None:
        yield header, "".join(seq_parts)


def already_done(out_tsv):
    done = set()
    if os.path.exists(out_tsv):
        with open(out_tsv) as f:
            next(f, None)  # header
            for line in f:
                mid = line.split("\t", 1)[0]
                if mid:
                    done.add(mid)
    return done


def parse_summary_line(stdout_text):
    # Last non-empty line is the one-line summary: "AGG x AMYLO y TURN z HELIX w HELAGG v BETA u"
    lines = [l for l in stdout_text.strip().splitlines() if l.strip()]
    if not lines:
        return None
    toks = lines[-1].split()
    vals = {}
    i = 0
    while i + 1 < len(toks):
        key, val = toks[i], toks[i + 1]
        if key in TANGO_FIELDS:
            vals[key] = val
        i += 2
    if set(vals) != set(TANGO_FIELDS):
        return None
    return vals


def main():
    chunk_fasta, out_tsv, tango_bin = sys.argv[1], sys.argv[2], sys.argv[3]
    work_dir = f"/tmp/tango_work_{uuid.uuid4().hex[:8]}"
    os.makedirs(work_dir, exist_ok=True)

    done = already_done(out_tsv)
    print(f"Resuming: {len(done)} sequences already done in {out_tsv}")

    new_file = not os.path.exists(out_tsv)
    out_f = open(out_tsv, "a")
    if new_file:
        out_f.write("mutation_id\t" + "\t".join(TANGO_FIELDS) + "\n")

    n_total = n_skipped = n_run = n_failed = 0
    for mutation_id, seq in read_fasta(chunk_fasta):
        n_total += 1
        if mutation_id in done:
            n_skipped += 1
            continue

        call_name = os.path.join(work_dir, "call")
        try:
            result = subprocess.run(
                [tango_bin, "call", "ct=N", "nt=N", "ph=7.4", "te=298", "io=0.05", f"seq={seq}"],
                cwd=work_dir,
                capture_output=True,
                text=True,
                timeout=1800,
            )
            vals = parse_summary_line(result.stdout)
        except subprocess.TimeoutExpired:
            vals = None
        finally:
            txt_path = call_name + ".txt"
            if os.path.exists(txt_path):
                os.remove(txt_path)

        if vals is None:
            n_failed += 1
            print(f"  FAILED to parse TANGO output for {mutation_id} (len={len(seq)})")
            continue

        out_f.write(mutation_id + "\t" + "\t".join(vals[k] for k in TANGO_FIELDS) + "\n")
        out_f.flush()
        n_run += 1
        if n_run % 5000 == 0:
            print(f"  {n_run} new sequences scored ({n_skipped} skipped, {n_failed} failed)")

    out_f.close()
    os.rmdir(work_dir) if os.path.isdir(work_dir) and not os.listdir(work_dir) else None

    print(f"\nTotal in chunk: {n_total}")
    print(f"Already done (skipped): {n_skipped}")
    print(f"Newly scored: {n_run}")
    print(f"Failed to parse: {n_failed}")


if __name__ == "__main__":
    main()
