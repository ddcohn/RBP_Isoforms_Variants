"""Step 2b: for RBPs whose scored transcript does NOT encode a protein of the input
length, measure how much of the input protein is still present in the scored protein.

`input_residues_in_scored_protein` = (residues of the input sequence that fall in matching
blocks against the scored protein, found with difflib.SequenceMatcher) / input length.
1.0 means every input residue is found in order in the scored protein (the scored isoform
can be longer); small values mean the scored isoform is a poor stand-in for the input.
It is a block-matching coverage, not an alignment identity.

Output: results/tables/isoform_similarity_mismatches.tsv
Run from the 13_evolutionary_conservation directory.
"""
import difflib
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import biomart as bm  # noqa: E402

INPUT = Path("../12_variant_effect_prediction/rbp_filter/rbp_wt_list_with_go.tsv")
OUT = Path("results/tables")


def main():
    m = pd.read_csv(OUT / "rbp_transcript_map.tsv", sep="\t")
    bad = m[(m["match_type"] != "none") & (~m["length_matches_input"])]
    seq = pd.read_csv(INPUT, sep="\t", dtype=str, usecols=["GeneSymbol", "UniProtID", "Sequence"])
    bad = bad.merge(seq, on=["GeneSymbol", "UniProtID"])
    pep = bm.query("hsapiens_gene_ensembl", ["ensembl_transcript_id", "peptide"],
                   {"ensembl_transcript_id": bad["ensembl_transcript_id"].tolist()})
    # BioMart returns the sequence column first; address columns by name
    pm = {r["Transcript stable ID"]: r["Peptide"].rstrip("*") for _, r in pep.iterrows()}
    rows = []
    for _, r in bad.iterrows():
        a, b = r["Sequence"], pm.get(r["ensembl_transcript_id"], "")
        sm = difflib.SequenceMatcher(None, a, b, autojunk=False)
        frac = sum(blk.size for blk in sm.get_matching_blocks()) / len(a)
        rows.append({"GeneSymbol": r["GeneSymbol"], "UniProtID": r["UniProtID"],
                     "input_aa": len(a), "scored_aa": len(b),
                     "input_residues_in_scored_protein": round(frac, 4)})
    out = pd.DataFrame(rows).sort_values("input_residues_in_scored_protein")
    out.to_csv(OUT / "isoform_similarity_mismatches.tsv", sep="\t", index=False)
    print(out.to_string(index=False))
    print(f"{(out['input_residues_in_scored_protein'] >= 0.9).sum()} of {len(out)} keep >= 90% of input residues")


if __name__ == "__main__":
    main()
