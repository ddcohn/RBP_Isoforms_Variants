"""Step 4: merge transcript mapping, CDS conservation and Compara orthology into
one row per RBP.  Run from the 13_evolutionary_conservation directory.

Output: results/tables/rbp_conservation_scores.tsv
Orthology columns (per species, from Ensembl Compara via BioMart):
  <sp>_orthology_type   best row per gene: one2one preferred, else highest query %id
  <sp>_pid_query        % of the human protein identical to the ortholog
  <sp>_pid_target       % of the ortholog identical to the human protein
  <sp>_goc, <sp>_wga    gene-order conservation score, whole-genome-alignment coverage
                        (only where Ensembl provides them)
  <sp>_n_orthologs      number of ortholog rows for the gene (1:many gives >1)
Mapping QC: length_matches_input, expanded_search, input_residues_in_scored_protein
  (1.0 when the scored protein has the input length; otherwise from step 2b).
Summary: n_species_with_ortholog, n_species_one2one, deepest_ortholog_species.
Absence of a Compara call is not proof the gene is absent from that species.
"""
from pathlib import Path

import numpy as np
import pandas as pd

RAW, OUT = Path("data/raw"), Path("results/tables")
INPUT = Path("../12_variant_effect_prediction/rbp_filter/rbp_wt_list_with_go.tsv")
# ordered close -> distant (divergence from human)
LADDER = ["ptroglodytes", "mmulatta", "mmusculus", "rnorvegicus", "ggallus",
          "xtropicalis", "drerio", "dmelanogaster", "celegans", "scerevisiae"]


def load_ortho(sp):
    df = pd.read_csv(RAW / f"ortho_{sp}.tsv", sep="\t", dtype=str)
    return df, df.columns.tolist()


def summarise(sp):
    df, hdr = load_ortho(sp)
    names = {}
    for h in hdr[1:]:
        hl = h.lower()
        if "stable id" in hl:
            names[h] = "gene"
        elif "homology type" in hl:
            names[h] = "orthology_type"
        elif "target" in hl and "identical to query" in hl:
            names[h] = "pid_target"
        elif "query" in hl and "identical to target" in hl:
            names[h] = "pid_query"
        elif "gene-order" in hl:
            names[h] = "goc"
        elif "whole-genome" in hl:
            names[h] = "wga"
        elif "confidence" in hl:
            names[h] = "confidence"
        else:
            raise ValueError(f"unrecognised BioMart column {h!r} for {sp}")
    df = df.rename(columns={hdr[0]: "ensembl_gene_id", **names})
    df = df.dropna(subset=["gene"])
    for c in ("pid_target", "pid_query", "goc", "wga", "confidence"):
        if c in df:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    df["is_o2o"] = df["orthology_type"].eq("ortholog_one2one")
    df = df.sort_values(["ensembl_gene_id", "is_o2o", "pid_query"],
                        ascending=[True, False, False])
    n = df.groupby("ensembl_gene_id").size().rename("n_orthologs")
    best = df.groupby("ensembl_gene_id", as_index=False).first().merge(n, on="ensembl_gene_id")
    keep = ["ensembl_gene_id", "orthology_type", "pid_query", "pid_target", "n_orthologs"]
    keep += [c for c in ("goc", "wga") if c in best]
    best = best[keep].rename(columns={c: f"{sp}_{c}" for c in keep if c != "ensembl_gene_id"})
    return best


def main():
    m = pd.read_csv(OUT / "rbp_transcript_map.tsv", sep="\t")
    cons = pd.read_csv(OUT / "cds_conservation_by_transcript.tsv", sep="\t")
    ann = pd.read_csv(INPUT, sep="\t", dtype=str,
                      usecols=["GeneSymbol", "UniProtID", "Ribosomal_protein",
                               "Spliceosomal_protein", "Other_RNP_machine", "HGNC_Gene_Family"])
    sim = pd.read_csv(OUT / "isoform_similarity_mismatches.tsv", sep="\t",
                      usecols=["GeneSymbol", "UniProtID", "input_residues_in_scored_protein"])
    m = m.merge(sim, on=["GeneSymbol", "UniProtID"], how="left")
    # exact-length matches contain the whole input protein by construction
    m.loc[m["length_matches_input"].eq(True), "input_residues_in_scored_protein"] = 1.0
    t = m.merge(cons, on="ensembl_transcript_id", how="left")
    t = t.merge(ann, on=["GeneSymbol", "UniProtID"], how="left")
    for sp in LADDER:
        t = t.merge(summarise(sp), on="ensembl_gene_id", how="left")
    has = pd.DataFrame({sp: t[f"{sp}_orthology_type"].notna() for sp in LADDER})
    o2o = pd.DataFrame({sp: t[f"{sp}_orthology_type"].eq("ortholog_one2one") for sp in LADDER})
    t["n_species_with_ortholog"] = has.sum(axis=1)
    t["n_species_one2one"] = o2o.sum(axis=1)
    deepest = []
    for _, row in has.iterrows():
        idx = [i for i, sp in enumerate(LADDER) if row[sp]]
        deepest.append(LADDER[max(idx)] if idx else "none")
    t["deepest_ortholog_species"] = deepest
    t.loc[t["match_type"].eq("none"), ["n_species_with_ortholog", "n_species_one2one",
                                      "deepest_ortholog_species"]] = [np.nan, np.nan, np.nan]
    t.to_csv(OUT / "rbp_conservation_scores.tsv", sep="\t", index=False)
    print(t.shape)
    print(t[["phyloP100way_mean", "phastCons100way_mean", "phyloP30way_mean",
             "phastCons30way_mean"]].describe().round(3).to_string())
    print("min CDS coverage (phyloP100way):", round(t["phyloP100way_cov"].min(), 3))


if __name__ == "__main__":
    main()
