"""Step 9: phyloP and phastCons for EVERY coding base of every selected transcript, plus a
per-residue summary.  Both tracks from both alignments (100 vertebrates, 30-way), read remotely
from UCSC hg38 bigWigs.

Outputs (results/per_base/, git-ignored because of size):
  rbp_cds_per_base_scores.parquet   one row per CDS base (about 5.0 million rows)
  rbp_per_residue_scores.tsv.gz     one row per amino-acid residue (codon means)

Bases are in transcript (5'->3') order.  For minus-strand genes the genomic positions of each exon
are reversed so cds_pos, codon_number and codon_position follow the reading frame.  All selected
transcripts have a CDS length divisible by 3, so the frame starts at base 1.
Missing alignment data stay empty (NaN); they are never zero-filled.
Run from the 13_evolutionary_conservation directory.
"""
import subprocess
import sys
import tempfile
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
T = Path("results/tables")
PB = Path("results/per_base")
_h = {}


def read(track, chrom, s0, e0, tries=6):
    for k in range(tries):
        try:
            if _h.get(track) is None:
                _h[track] = pyBigWig.open(TRACKS[track])
            return np.asarray(_h[track].values(chrom, s0, e0), dtype=np.float32)
        except Exception:
            _h[track] = None
            time.sleep(5 * (k + 1))
    raise RuntimeError(f"failed {track} {chrom}:{s0}-{e0}")


def one_transcript(tx, segs, gene, uniprot, peptide):
    segs = segs.sort_values("exon_rank")
    chrom, strand = segs["chrom"].iloc[0], int(segs["strand"].iloc[0])
    lo = int(segs["cds_start_1based"].min()) - 1
    hi = int(segs["cds_end_1based"].max())
    span = {t: read(t, chrom, lo, hi) for t in TRACKS}
    pos, exon, vals = [], [], {t: [] for t in TRACKS}
    for r in segs.itertuples():
        p = np.arange(r.cds_start_1based, r.cds_end_1based + 1)
        sl = slice(r.cds_start_1based - 1 - lo, r.cds_end_1based - lo)
        v = {t: span[t][sl] for t in TRACKS}
        if strand == -1:
            p, v = p[::-1], {t: x[::-1] for t, x in v.items()}
        pos.append(p)
        exon.append(np.full(len(p), r.exon_rank))
        for t in TRACKS:
            vals[t].append(v[t])
    pos, exon = np.concatenate(pos), np.concatenate(exon)
    n = len(pos)
    assert n % 3 == 0 and n // 3 - 1 == len(peptide), (tx, n, len(peptide))
    cds_pos = np.arange(1, n + 1)
    codon = (cds_pos - 1) // 3 + 1
    cpos = ((cds_pos - 1) % 3 + 1).astype(np.int8)
    aa = np.array(list(peptide) + ["*"])[codon - 1]
    df = pd.DataFrame({"GeneSymbol": gene, "ensembl_transcript_id": tx, "chrom": chrom,
                       "strand": np.int8(strand), "genomic_pos": pos, "exon_rank": exon.astype(np.int16),
                       "cds_pos": cds_pos, "codon_number": codon, "codon_position": cpos, "aa": aa})
    for t in TRACKS:
        df[t] = np.concatenate(vals[t])
    # residue summary (stop codon excluded)
    res = df[df["aa"] != "*"]
    g = res.groupby("codon_number")
    out = pd.DataFrame({"GeneSymbol": gene, "UniProtID": uniprot, "ensembl_transcript_id": tx,
                        "residue_number": g.size().index.to_numpy(),
                        "aa": g["aa"].first().to_numpy()})
    p12 = res[res["codon_position"] < 3].groupby("codon_number")
    for t in TRACKS:
        out[f"{t}_codon_mean"] = g[t].mean().to_numpy().round(3)
        out[f"{t}_pos12_mean"] = p12[t].mean().reindex(g.size().index).to_numpy().round(3)
    return df, out


def load_inputs():
    m = pd.read_csv(T / "rbp_transcript_map.tsv", sep="\t")
    m = m[m["match_type"] != "none"].sort_values("ensembl_transcript_id")
    seg = pd.read_csv(T / "rbp_cds_segments.tsv", sep="\t")
    pep = pd.read_csv("data/raw/peptides.tsv", sep="\t", dtype=str).set_index("tx")["peptide"]
    return m, seg, pep


def run_shard(i, n, tmp):
    m, seg, pep = load_inputs()
    seg_by_tx = {k: g for k, g in seg.groupby("ensembl_transcript_id")}
    bases, resid = [], []
    for k, r in enumerate(m.iloc[i::n].itertuples(), 1):
        b, rs = one_transcript(r.ensembl_transcript_id, seg_by_tx[r.ensembl_transcript_id],
                               r.GeneSymbol, r.UniProtID, pep[r.ensembl_transcript_id])
        bases.append(b)
        resid.append(rs)
        if k % 100 == 0:
            print(f"shard {i}: {k}", flush=True)
    pd.concat(bases, ignore_index=True).to_parquet(Path(tmp) / f"b_{i}.parquet", compression="zstd")
    pd.concat(resid, ignore_index=True).to_parquet(Path(tmp) / f"r_{i}.parquet", compression="zstd")


def main(workers=4):
    PB.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    with tempfile.TemporaryDirectory() as tmp:
        procs = [subprocess.Popen([sys.executable, __file__, "--shard", str(i), str(workers), tmp])
                 for i in range(workers)]
        codes = [p.wait() for p in procs]
        assert all(c == 0 for c in codes), f"worker exit codes {codes}"
        base = pd.concat([pd.read_parquet(Path(tmp) / f"b_{i}.parquet") for i in range(workers)],
                         ignore_index=True)
        res = pd.concat([pd.read_parquet(Path(tmp) / f"r_{i}.parquet") for i in range(workers)],
                        ignore_index=True)
    base = base.sort_values(["GeneSymbol", "ensembl_transcript_id", "cds_pos"], ignore_index=True)
    res = res.sort_values(["GeneSymbol", "ensembl_transcript_id", "residue_number"], ignore_index=True)

    # --- checks -------------------------------------------------------------------------
    seg = pd.read_csv(T / "rbp_cds_segments.tsv", sep="\t")
    assert len(base) == seg["cds_len"].sum(), "row count != total CDS bases"
    assert not base.duplicated(["ensembl_transcript_id", "cds_pos"]).any()
    old = pd.read_csv(T / "cds_conservation_by_transcript.tsv", sep="\t").set_index("ensembl_transcript_id")
    for t in TRACKS:
        new = base.groupby("ensembl_transcript_id")[t].mean()
        d = (new - old[f"{t}_mean"].reindex(new.index)).abs().max()
        print(f"max |per-base mean - step-3 mean| {t}: {d:.2e}", flush=True)
        assert d < 1e-3
    chk = base.groupby(["strand", "codon_position"])["phyloP100way"].mean().round(2).unstack()
    print("mean phyloP100way by strand x codon position:\n", chk.to_string(), flush=True)
    assert (chk[3] < chk[1]).all() and (chk[3] < chk[2]).all(), "codon position 3 is not lowest: frame/strand error"

    base.to_parquet(PB / "rbp_cds_per_base_scores.parquet", compression="zstd", index=False)
    res.to_csv(PB / "rbp_per_residue_scores.tsv.gz", sep="\t", index=False, compression="gzip")
    print(f"per-base rows {len(base):,}; per-residue rows {len(res):,}; {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    if len(sys.argv) == 5 and sys.argv[1] == "--shard":
        run_shard(int(sys.argv[2]), int(sys.argv[3]), sys.argv[4])
    else:
        main()
