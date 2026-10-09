"""Step 7: one row per RBP, one headline conservation number, plus a reliability flag.

Output: results/tables/rbp_conservation_simple.tsv
Headline score = mean phyloP over the CDS, 100-vertebrate alignment (`cds_phyloP_100way`).
Percentile and tier are relative to the mapped RBPs in this list, not absolute thresholds.
Run from the 13_evolutionary_conservation directory.
"""
from pathlib import Path

import numpy as np
import pandas as pd

T = Path("results/tables")


def main():
    t = pd.read_csv(T / "rbp_conservation_scores.tsv", sep="\t")
    mapped = t["match_type"].ne("none")
    s = pd.DataFrame({"GeneSymbol": t["GeneSymbol"], "UniProtID": t["UniProtID"]})
    s["cds_phyloP_100way"] = t["phyloP100way_mean"].round(3)
    s["conservation_percentile"] = (t["phyloP100way_mean"].where(mapped).rank(pct=True) * 100).round(1)
    s["conservation_tier"] = np.select(
        [s["conservation_percentile"] >= 75, s["conservation_percentile"] <= 25],
        ["high", "low"], default="medium")
    s.loc[~mapped, "conservation_tier"] = ""
    s["mouse_protein_identity_pct"] = t["mmusculus_pid_query"].round(1)
    s["n_species_with_ortholog"] = t["n_species_with_ortholog"]
    s["deepest_ortholog_species"] = t["deepest_ortholog_species"]

    iso_ok = t["input_residues_in_scored_protein"] >= 0.9
    cov_ok = t["phyloP100way_cov"] >= 0.9
    s["score_reliable"] = mapped & iso_ok & cov_ok
    notes = []
    for m, i, c in zip(mapped, iso_ok, cov_ok):
        n = []
        if not m:
            n.append("not mapped to an Ensembl transcript")
        else:
            if not i:
                n.append("scored isoform covers <90% of the input protein")
            if not c:
                n.append("<90% of CDS bases have alignment data")
        notes.append("; ".join(n))
    s["reliability_note"] = notes
    s.to_csv(T / "rbp_conservation_simple.tsv", sep="\t", index=False)
    print(s["conservation_tier"].value_counts().to_dict(), "| reliable:", int(s["score_reliable"].sum()),
          "of", len(s))
    print(s.loc[s["GeneSymbol"].isin(["HNRNPA1", "SRSF1", "ELAVL1", "NXF2", "UNKL"])].to_string(index=False))


if __name__ == "__main__":
    main()
