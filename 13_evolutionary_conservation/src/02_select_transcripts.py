"""Step 2: map each RBP (GeneSymbol + UniProtID) to one Ensembl transcript and
fetch that transcript's CDS exon segments.

Transcript choice per row:
  1. candidate transcripts = protein-coding transcripts annotated with the row's
     UniProt accession (match_type 'uniprot'); if none, transcripts of genes whose
     symbol equals GeneSymbol (match_type 'symbol');
  1b. if no candidate encodes a protein of exactly the input `Length`, widen the search to
     every protein-coding transcript of the candidates' gene(s) plus the Ensembl transcripts
     cross-referenced from the UniProt entry, and use those only if one matches exactly
     (column `expanded_search`);
  2. among candidates prefer (a) a transcript whose CDS encodes a protein of exactly the
     input `Length` (so the scored CDS is the same protein as the input sequence), then
     (b) MANE Select, (c) Ensembl canonical, (d) longest CDS, (e) lowest ENST.

Outputs:
  results/tables/rbp_transcript_map.tsv   one row per input RBP
  results/tables/rbp_cds_segments.tsv     one row per CDS exon segment of the chosen transcript
Run from the 13_evolutionary_conservation directory.
"""
import sys
from pathlib import Path

import pandas as pd
import requests

sys.path.insert(0, str(Path(__file__).parent))
import biomart as bm  # noqa: E402

INPUT = Path("../12_variant_effect_prediction/rbp_filter/rbp_wt_list_with_go.tsv")
RAW = Path("data/raw")
OUT = Path("results/tables")
OUT.mkdir(parents=True, exist_ok=True)
STD_CHROMS = {str(i) for i in range(1, 23)} | {"X", "Y", "MT"}


def load_tx_table():
    nu = pd.read_csv(RAW / "tx_map_names_uniprot.tsv", sep="\t", dtype=str)
    nu.columns = ["gene_id", "tx_id", "symbol", "uniprot", "chrom"]
    can = pd.read_csv(RAW / "tx_map_canonical.tsv", sep="\t", dtype=str)
    can.columns = ["gene_id", "tx_id", "is_canonical"]
    mane = pd.read_csv(RAW / "tx_map_mane.tsv", sep="\t", dtype=str)
    mane.columns = ["gene_id", "tx_id", "mane_select"]
    cds = pd.read_csv(RAW / "tx_map_cds_length.tsv", sep="\t", dtype=str)
    cds.columns = ["gene_id", "tx_id", "cds_length"]
    base = nu[["gene_id", "tx_id", "symbol", "chrom"]].drop_duplicates()
    # symbol/chrom are constant per transcript; uniprot is the multi-valued column
    base = base.groupby(["gene_id", "tx_id"], as_index=False).first()
    for extra in (can.drop_duplicates(), mane.drop_duplicates(), cds.drop_duplicates()):
        base = base.merge(extra, on=["gene_id", "tx_id"], how="left")
    base["is_canonical"] = base["is_canonical"].eq("1")
    base["is_mane"] = base["mane_select"].notna()
    base["cds_length"] = pd.to_numeric(base["cds_length"], errors="coerce")
    base = base[base["chrom"].isin(STD_CHROMS) & base["cds_length"].notna()]
    uni = nu.dropna(subset=["uniprot"])[["tx_id", "uniprot"]].drop_duplicates()
    return base, uni


def rank(df, input_len):
    df = df.assign(len_match=(df["cds_length"] / 3 - 1) == float(input_len))
    return df.sort_values(
        ["len_match", "is_mane", "is_canonical", "cds_length", "tx_id"],
        ascending=[False, False, False, False, True],
    )


def uniprot_ensembl_xrefs(accessions):
    """UniProt accession -> set of Ensembl transcript IDs (version stripped) listed in the entry."""
    import io
    import re
    out = {}
    accs = sorted(set(accessions))
    for i in range(0, len(accs), 200):
        r = requests.get("https://rest.uniprot.org/uniprotkb/accessions",
                         params={"accessions": ",".join(accs[i:i + 200]),
                                 "fields": "accession,xref_ensembl", "format": "tsv", "size": 500},
                         timeout=120)
        r.raise_for_status()
        df = pd.read_csv(io.StringIO(r.text), sep="\t", dtype=str).fillna("")
        for _, row in df.iterrows():
            out[row["Entry"]] = set(re.findall(r"(ENST\d+)\.\d+", row["Ensembl"]))
    return out


def main():
    rbp = pd.read_csv(INPUT, sep="\t", dtype=str, usecols=["GeneSymbol", "UniProtID", "Length"])
    tx, uni = load_tx_table()
    by_uni = uni.merge(tx, on="tx_id")
    by_sym = tx.copy()

    # first pass: which rows have no candidate transcript of exactly the input length?
    def first_cands(r):
        c, h = by_uni[by_uni["uniprot"] == r["UniProtID"]], "uniprot"
        if c.empty:
            c, h = by_sym[by_sym["symbol"] == r["GeneSymbol"]], "symbol"
        return c, h

    need = [r["UniProtID"] for _, r in rbp.iterrows()
            if not ((lambda c: (c["cds_length"] / 3 - 1 == float(r["Length"])).any())(first_cands(r)[0]))]
    xrefs = uniprot_ensembl_xrefs(need)
    print(f"{len(need)} rows lack an exact-length candidate; widening search for them", flush=True)

    rows = []
    for _, r in rbp.iterrows():
        cand, how = first_cands(r)
        if cand.empty:
            rows.append({"GeneSymbol": r["GeneSymbol"], "UniProtID": r["UniProtID"],
                         "input_length": r["Length"], "match_type": "none"})
            continue
        expanded = False
        if not (cand["cds_length"] / 3 - 1 == float(r["Length"])).any():
            ext = tx[(tx["gene_id"].isin(set(cand["gene_id"])) | tx["tx_id"].isin(xrefs.get(r["UniProtID"], set())))
                     & (tx["cds_length"] / 3 - 1 == float(r["Length"]))]
            if not ext.empty:
                cand, expanded = ext, True
        best = rank(cand, r["Length"]).iloc[0]
        rule = "MANE_Select" if best["is_mane"] else "Ensembl_canonical" if best["is_canonical"] else "longest_CDS"
        rows.append({
            "GeneSymbol": r["GeneSymbol"], "UniProtID": r["UniProtID"],
            "input_length": r["Length"], "match_type": how,
            "ensembl_gene_id": best["gene_id"], "ensembl_transcript_id": best["tx_id"],
            "ensembl_symbol": best["symbol"], "chrom": "chr" + best["chrom"],
            "cds_length_nt": int(best["cds_length"]), "selection_rule": rule,
            "n_candidate_transcripts": cand["tx_id"].nunique(),
            "expanded_search": expanded,
        })
    m = pd.DataFrame(rows)
    m["symbol_agrees"] = m["GeneSymbol"] == m["ensembl_symbol"]
    # QC against the user's own protein length: CDS/3 - 1 stop codon
    m["cds_aa"] = m["cds_length_nt"] / 3 - 1
    m["length_matches_input"] = m["cds_aa"] == pd.to_numeric(m["input_length"])
    m.to_csv(OUT / "rbp_transcript_map.tsv", sep="\t", index=False)
    print(m["match_type"].value_counts().to_dict(), flush=True)
    print(m["selection_rule"].value_counts().to_dict(), flush=True)
    print("symbol agrees:", int(m["symbol_agrees"].sum()),
          "| CDS-derived length == input length:", int(m["length_matches_input"].sum()), flush=True)

    # CDS segments for the chosen transcripts, in chunks (GET URL length limit)
    txs = sorted(m["ensembl_transcript_id"].dropna().unique())
    attrs = ["ensembl_transcript_id", "ensembl_exon_id", "rank", "chromosome_name",
             "strand", "genomic_coding_start", "genomic_coding_end"]
    parts = []
    for i in range(0, len(txs), 100):
        parts.append(bm.query("hsapiens_gene_ensembl", attrs,
                              {"ensembl_transcript_id": txs[i:i + 100]}))
        print(f"structure chunk {i // 100 + 1}/{(len(txs) + 99) // 100}", flush=True)
    seg = pd.concat(parts, ignore_index=True)
    seg.columns = ["ensembl_transcript_id", "exon_id", "exon_rank", "chrom", "strand",
                   "cds_start_1based", "cds_end_1based"]
    seg = seg.dropna(subset=["cds_start_1based", "cds_end_1based"]).copy()
    for c in ("exon_rank", "strand", "cds_start_1based", "cds_end_1based"):
        seg[c] = pd.to_numeric(seg[c])
    seg["chrom"] = "chr" + seg["chrom"].astype(str)
    seg["cds_len"] = seg["cds_end_1based"] - seg["cds_start_1based"] + 1
    seg = seg.sort_values(["ensembl_transcript_id", "exon_rank"])
    seg.to_csv(OUT / "rbp_cds_segments.tsv", sep="\t", index=False)
    tot = seg.groupby("ensembl_transcript_id")["cds_len"].sum()
    chk = m.dropna(subset=["ensembl_transcript_id"]).set_index("ensembl_transcript_id")["cds_length_nt"]
    agree = (tot.reindex(chk.index) == chk).sum()
    print(f"segments: {len(seg)}; transcripts whose summed CDS segments == BioMart cds_length: {agree}/{len(chk)}", flush=True)


if __name__ == "__main__":
    main()
