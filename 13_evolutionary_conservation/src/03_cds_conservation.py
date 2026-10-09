"""Step 3: mean CDS conservation per transcript and per CDS exon segment from
UCSC hg38 phyloP / phastCons bigWigs, read remotely (HTTP range requests; no
multi-GB download).

Tracks (UCSC goldenPath/hg38):
  phyloP100way / phastCons100way : 100-way vertebrate alignment
  phyloP30way  / phastCons30way  : 30-way alignment (27 primates + 3 outgroups)

Per transcript: pooled over every CDS base, so exons are length-weighted.
Bases with no alignment data are excluded from the mean and reported through
`cov_<track>` (covered CDS bases / total CDS bases) -- they are NOT zero-filled.
Run from the 13_evolutionary_conservation directory.
"""
import time
from pathlib import Path

import numpy as np
import pandas as pd
import pyBigWig

BASE = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38"
TRACKS = {
    "phyloP100way": f"{BASE}/phyloP100way/hg38.phyloP100way.bw",
    "phastCons100way": f"{BASE}/phastCons100way/hg38.phastCons100way.bw",
    "phyloP30way": f"{BASE}/phyloP30way/hg38.phyloP30way.bw",
    "phastCons30way": f"{BASE}/phastCons30way/hg38.phastCons30way.bw",
}
# "conserved base" cutoffs for the fraction column
CUTOFF = {"phyloP": 2.0, "phastCons": 0.9}
OUT = Path("results/tables")
_handles = {}


def handle(track):
    if track not in _handles or _handles[track] is None:
        _handles[track] = pyBigWig.open(TRACKS[track])
    return _handles[track]


def read(track, chrom, s0, e0, tries=6):
    for k in range(tries):
        try:
            v = handle(track).values(chrom, s0, e0)
            return np.asarray(v, dtype=float)
        except Exception:
            _handles[track] = None
            time.sleep(5 * (k + 1))
    raise RuntimeError(f"failed {track} {chrom}:{s0}-{e0}")


def work(item):
    tx, segs = item
    chrom = segs[0][2]
    lo = min(s1 for (_, _, _, s1, _) in segs) - 1
    hi = max(e1 for (_, _, _, _, e1) in segs)
    seg_rows, tx_row = [], {"ensembl_transcript_id": tx}
    pooled = {t: [] for t in TRACKS}
    # one request per track spanning the whole CDS extent (far fewer HTTP range
    # requests than one per exon); CDS segments are then sliced out locally
    span = {t: read(t, chrom, lo, hi) for t in TRACKS}
    for (exon_id, rank, chrom, s1, e1) in segs:
        s0, e0 = s1 - 1, e1  # 1-based inclusive -> 0-based half-open
        row = {"ensembl_transcript_id": tx, "exon_id": exon_id, "exon_rank": rank,
               "chrom": chrom, "start0": s0, "end": e0, "n_bases": e0 - s0}
        for t in TRACKS:
            v = span[t][s0 - lo:e0 - lo]
            pooled[t].append(v)
            ok = v[np.isfinite(v)]
            row[f"{t}_mean"] = float(ok.mean()) if ok.size else np.nan
            row[f"{t}_covered"] = int(ok.size)
        seg_rows.append(row)
    for t, arrs in pooled.items():
        v = np.concatenate(arrs)
        ok = v[np.isfinite(v)]
        fam = "phyloP" if t.startswith("phyloP") else "phastCons"
        tx_row[f"{t}_mean"] = float(ok.mean()) if ok.size else np.nan
        tx_row[f"{t}_median"] = float(np.median(ok)) if ok.size else np.nan
        tx_row[f"{t}_frac_conserved"] = float((ok >= CUTOFF[fam]).mean()) if ok.size else np.nan
        tx_row[f"{t}_cov"] = float(ok.size / v.size) if v.size else np.nan
    tx_row["cds_bases"] = int(sum(len(a) for a in pooled["phyloP100way"]))
    return tx_row, seg_rows


def load_items():
    seg = pd.read_csv(OUT / "rbp_cds_segments.tsv", sep="\t")
    items = []
    for tx, g in seg.groupby("ensembl_transcript_id", sort=True):
        items.append((tx, [(r.exon_id, int(r.exon_rank), r.chrom, int(r.cds_start_1based),
                            int(r.cds_end_1based)) for r in g.itertuples()]))
    return items


def run_shard(i, n):
    """Worker: process every n-th transcript starting at i; write a shard file."""
    items = load_items()[i::n]
    tx_rows, seg_rows = [], []
    for k, item in enumerate(items, 1):
        trow, srows = work(item)
        tx_rows.append(trow)
        seg_rows.extend(srows)
        if k % 100 == 0:
            print(f"shard {i}: {k}/{len(items)}", flush=True)
    SH = OUT / "shards"
    SH.mkdir(exist_ok=True)
    pd.DataFrame(tx_rows).to_csv(SH / f"tx_{i}.tsv", sep="\t", index=False)
    pd.DataFrame(seg_rows).to_csv(SH / f"seg_{i}.tsv", sep="\t", index=False)


def main(workers=4):
    # Independent subprocesses (the sandbox forbids multiprocessing semaphores).
    import subprocess
    import sys
    t0 = time.time()
    procs = [subprocess.Popen([sys.executable, __file__, "--shard", str(i), str(workers)])
             for i in range(workers)]
    codes = [p.wait() for p in procs]
    assert all(c == 0 for c in codes), f"worker exit codes {codes}"
    SH = OUT / "shards"
    tx = pd.concat([pd.read_csv(SH / f"tx_{i}.tsv", sep="\t") for i in range(workers)])
    sg = pd.concat([pd.read_csv(SH / f"seg_{i}.tsv", sep="\t") for i in range(workers)])
    tx.sort_values("ensembl_transcript_id").to_csv(OUT / "cds_conservation_by_transcript.tsv", sep="\t", index=False)
    sg.sort_values(["ensembl_transcript_id", "exon_rank"]).to_csv(OUT / "cds_conservation_by_exon.tsv", sep="\t", index=False)
    print(f"done {len(tx)} transcripts, {len(sg)} segments in {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    import sys
    if len(sys.argv) == 4 and sys.argv[1] == "--shard":
        run_shard(int(sys.argv[2]), int(sys.argv[3]))
    else:
        main()
