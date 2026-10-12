"""Step 8: fetch the Ensembl protein sequence of every selected transcript, so each coding
codon can be labelled with its amino acid.

Output: data/raw/peptides.tsv (transcript, peptide).  Also checks each peptide's length against
the transcript's CDS length and, for exact-length matches, how often the peptide equals the input
list's UniProt sequence.  Run from the 13_evolutionary_conservation directory.
"""
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import biomart as bm  # noqa: E402

INPUT = Path("../12_variant_effect_prediction/rbp_filter/rbp_wt_list_with_go.tsv")
OUT = Path("data/raw/peptides.tsv")


def main():
    m = pd.read_csv("results/tables/rbp_transcript_map.tsv", sep="\t")
    m = m[m["match_type"] != "none"]
    txs = sorted(m["ensembl_transcript_id"].unique())
    parts = []
    for i in range(0, len(txs), 150):
        p = bm.query("hsapiens_gene_ensembl", ["ensembl_transcript_id", "peptide"],
                     {"ensembl_transcript_id": txs[i:i + 150]})
        parts.append(p.rename(columns={"Transcript stable ID": "tx", "Peptide": "peptide"}))
        print(f"peptide chunk {i // 150 + 1}/{(len(txs) + 149) // 150}", flush=True)
    pep = pd.concat(parts, ignore_index=True)
    pep["peptide"] = pep["peptide"].str.rstrip("*")
    pep = pep.drop_duplicates("tx")
    pep.to_csv(OUT, sep="\t", index=False)

    chk = m.merge(pep, left_on="ensembl_transcript_id", right_on="tx", how="left")
    assert chk["peptide"].notna().all(), "missing peptides: " + str(chk.loc[chk["peptide"].isna(), "GeneSymbol"].tolist()[:10])
    assert (chk["peptide"].str.len() == chk["cds_aa"]).all(), "peptide length != CDS length / 3 - 1"
    seq = pd.read_csv(INPUT, sep="\t", dtype=str, usecols=["GeneSymbol", "UniProtID", "Sequence"])
    chk = chk.merge(seq, on=["GeneSymbol", "UniProtID"])
    ex = chk[chk["length_matches_input"]]
    same = (ex["peptide"] == ex["Sequence"]).sum()
    print(f"{len(chk)} peptides; lengths all equal CDS/3-1; exact-length rows whose peptide equals the "
          f"input sequence: {same}/{len(ex)}", flush=True)


if __name__ == "__main__":
    main()
