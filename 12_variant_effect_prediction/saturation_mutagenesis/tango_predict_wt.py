"""
Parallel version of the WT TANGO run: unlike the satmut chunk script
(which processes 140,000 sequences per single-core task, with the
array dimension providing the only parallelism), this is only 2,335
sequences total -- small enough that a ProcessPoolExecutor across
multiple cores on one node is the simpler way to go fast, rather than
splitting into an array job.

Each TANGO call is already independent (its own subprocess, own work
dir) so this parallelizes cleanly -- no shared state between workers
except writing finished results back to the main process to serialize
to disk.
"""
import csv
import os
import subprocess
import uuid
from concurrent.futures import ProcessPoolExecutor, as_completed

csv.field_size_limit(2**31 - 1)

TANGO_FIELDS = ["AGG", "AMYLO", "TURN", "HELIX", "HELAGG", "BETA"]
WT_LIST = "/u/project/kappel/ddcohn/protein_variant_effects/rbp_only/rbp_wt_list_with_go.tsv"
TANGO_BIN = "/u/home/d/ddcohn/tango_x86_64_release"
OUT_TSV = "/u/project/kappel/ddcohn/protein_variant_effects/rbp_only/wt_tango_scores.tsv"
N_WORKERS = 8


def parse_summary_line(stdout_text):
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


def score_one(acc, seq):
    work_dir = f"/tmp/tango_wt_work_{uuid.uuid4().hex[:8]}"
    os.makedirs(work_dir, exist_ok=True)
    call_name = os.path.join(work_dir, "call")
    try:
        result = subprocess.run(
            [TANGO_BIN, "call", "ct=N", "nt=N", "ph=7.4", "te=298", "io=0.05", f"seq={seq}"],
            cwd=work_dir,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True,
            timeout=1800,
        )
        vals = parse_summary_line(result.stdout)
    except subprocess.TimeoutExpired:
        vals = None
    finally:
        txt_path = call_name + ".txt"
        if os.path.exists(txt_path):
            os.remove(txt_path)
        if os.path.isdir(work_dir) and not os.listdir(work_dir):
            os.rmdir(work_dir)
    return acc, vals


def main():
    done = set()
    if os.path.exists(OUT_TSV):
        with open(OUT_TSV) as f:
            next(f, None)
            for line in f:
                a = line.split("\t", 1)[0]
                if a:
                    done.add(a)
    print(f"Resuming: {len(done)} WT proteins already done")

    todo = []
    with open(WT_LIST, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            acc = row["UniProtID"]
            if acc not in done:
                todo.append((acc, row["Sequence"]))

    print(f"To score: {len(todo)} proteins across {N_WORKERS} workers")

    new_file = not os.path.exists(OUT_TSV)
    out_f = open(OUT_TSV, "a")
    if new_file:
        out_f.write("uniprot_accession\t" + "\t".join(TANGO_FIELDS) + "\n")

    n_run = n_failed = 0
    with ProcessPoolExecutor(max_workers=N_WORKERS) as ex:
        futures = {ex.submit(score_one, acc, seq): acc for acc, seq in todo}
        for fut in as_completed(futures):
            acc, vals = fut.result()
            if vals is None:
                n_failed += 1
                print(f"  FAILED to parse TANGO output for {acc}")
                continue
            out_f.write(acc + "\t" + "\t".join(vals[k] for k in TANGO_FIELDS) + "\n")
            out_f.flush()
            n_run += 1
            if n_run % 200 == 0:
                print(f"  {n_run}/{len(todo)} scored ({n_failed} failed)")

    out_f.close()
    print(f"\nNewly scored: {n_run}")
    print(f"Failed to parse: {n_failed}")


if __name__ == "__main__":
    main()
