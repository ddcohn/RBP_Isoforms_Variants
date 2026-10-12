"""Step 6: write COLUMNS.md and results/tables/column_dictionary.tsv, a description of
every column in every results table.  The script fails if any column in a table lacks a
description (or a description names a column that no longer exists), so the dictionary
cannot silently drift from the tables.  Run from the 13_evolutionary_conservation directory.
"""
from pathlib import Path

import pandas as pd

T = Path("results/tables")
TRACKS = {
    "phyloP100way": "phyloP (100 vertebrates)",
    "phastCons100way": "phastCons (100 vertebrates)",
    "phyloP30way": "phyloP (30-way: 27 primates + 3 outgroups)",
    "phastCons30way": "phastCons (30-way: 27 primates + 3 outgroups)",
}
SPECIES = {
    "ptroglodytes": "chimpanzee", "mmulatta": "rhesus macaque", "mmusculus": "mouse",
    "rnorvegicus": "rat", "ggallus": "chicken", "xtropicalis": "Xenopus tropicalis (frog)",
    "drerio": "zebrafish", "dmelanogaster": "fruit fly", "celegans": "C. elegans (worm)",
    "scerevisiae": "budding yeast",
}
HAS_GOC = {"ptroglodytes", "mmulatta", "mmusculus", "rnorvegicus", "ggallus", "xtropicalis", "drerio"}
HAS_WGA = HAS_GOC - {"xtropicalis"}

TRACK_NOTE = {
    "phyloP": "phyloP is a per-base score on a -log10 p-value scale: positive = conserved "
              "(evolving more slowly than neutral), negative = accelerated, 0 = neutral",
    "phastCons": "phastCons is a per-base conservation probability from 0 to 1",
}


def track_cols(prefix=""):
    """Per-track columns shared by the by-transcript table and the main table."""
    out = []
    for tr, label in TRACKS.items():
        fam = "phyloP" if tr.startswith("phyloP") else "phastCons"
        cutoff = "phyloP >= 2" if fam == "phyloP" else "phastCons >= 0.9"
        out += [
            (f"{tr}_mean", f"Mean {label} over all CDS bases of the selected transcript "
                           f"(all exons pooled, so exons are length-weighted; bases without alignment "
                           f"data are excluded, not zero-filled). {TRACK_NOTE[fam]}."),
            (f"{tr}_median", f"Median {label} over the same CDS bases."),
            (f"{tr}_frac_conserved", f"Fraction of covered CDS bases with {cutoff} (0-1)."),
            (f"{tr}_cov", f"Fraction of all CDS bases that have {label} data (0-1); "
                          f"low values mean the mean/median rest on fewer bases."),
        ]
    return out


MAP_COLS = [
    ("GeneSymbol", "Gene symbol as given in the input RBP list."),
    ("UniProtID", "UniProt accession as given in the input RBP list."),
    ("input_length", "Protein length in amino acids from the input list (`Length`, the UniProt sequence)."),
    ("match_type", "How the row was linked to Ensembl: `uniprot` = a protein-coding transcript carries "
                   "this UniProt accession in Ensembl BioMart; `symbol` = no such transcript, so linked by "
                   "gene symbol; `none` = not linked (every score column is then empty)."),
    ("ensembl_gene_id", "Ensembl gene ID (release 116) of the selected transcript."),
    ("ensembl_transcript_id", "Ensembl transcript ID whose CDS was scored (version suffix omitted)."),
    ("ensembl_symbol", "Gene name Ensembl uses for the selected gene (may differ from `GeneSymbol` "
                       "for renamed genes)."),
    ("chrom", "GRCh38 chromosome of the selected transcript, `chr`-prefixed."),
    ("cds_length_nt", "Coding-sequence length of the selected transcript in nucleotides, including the "
                      "stop codon (from BioMart)."),
    ("selection_rule", "Which flag the selected transcript carries: `MANE_Select`, `Ensembl_canonical`, or "
                       "`longest_CDS` (neither flag; longest CDS among the candidates). An exact protein-length "
                       "match to the input is applied first, so this describes the transcript chosen, not "
                       "why it won."),
    ("n_candidate_transcripts", "Number of transcripts considered for this RBP (after any widening, see "
                                "`expanded_search`)."),
    ("expanded_search", "True if no initial candidate had the input protein length and the search was widened "
                        "to all coding transcripts of the gene plus the UniProt-cross-referenced Ensembl "
                        "transcripts, which produced an exact-length match. False otherwise."),
    ("symbol_agrees", "True if `GeneSymbol` equals `ensembl_symbol`."),
    ("cds_aa", "Length in amino acids of the protein encoded by the selected transcript "
               "(`cds_length_nt / 3 - 1`)."),
    ("length_matches_input", "True if `cds_aa` equals `input_length`, i.e. the scored CDS encodes a protein "
                             "of the same length as the input sequence. False for unmapped rows and for 28 "
                             "mapped RBPs with no exact-length transcript in Ensembl 116."),
]

SIM_COL = ("input_residues_in_scored_protein",
           "Fraction of the input protein's residues found, in order, in the scored protein (difflib "
           "matching blocks; a coverage measure, not an alignment identity). 1.0 for every exact-length "
           "match; for the 28 mismatches it ranges from 0.38 to 1.0 (1.0 means the scored isoform "
           "contains the whole input sequence and is longer). Empty for unmapped rows.")

ANNOT_COLS = [
    ("Ribosomal_protein", "`Y`/`N` flag carried over from the input list: ribosomal protein."),
    ("Spliceosomal_protein", "`Y`/`N` flag carried over from the input list: spliceosomal protein."),
    ("Other_RNP_machine", "`Y`/`N` flag carried over from the input list: member of another named RNP machine."),
    ("HGNC_Gene_Family", "HGNC gene-family name carried over from the input list (empty if none)."),
]


def ortho_cols():
    out = []
    for sp, name in SPECIES.items():
        out += [
            (f"{sp}_orthology_type", f"Compara homology type to {name} for the best row of this gene "
                                     f"(`ortholog_one2one`, `ortholog_one2many`, `ortholog_many2many`); "
                                     f"empty = no ortholog call. One-to-one is preferred, otherwise the row "
                                     f"with the highest `{sp}_pid_query`."),
            (f"{sp}_pid_query", f"Percent of the human protein that is identical to the {name} ortholog "
                                f"(Ensembl `%id. query gene identical to target gene`)."),
            (f"{sp}_pid_target", f"Percent of the {name} ortholog that is identical to the human protein "
                                 f"(Ensembl `%id. target gene identical to query gene`)."),
            (f"{sp}_n_orthologs", f"Number of {name} ortholog rows Ensembl lists for this human gene "
                                  f"(>1 for one-to-many and many-to-many)."),
        ]
        if sp in HAS_GOC:
            out.append((f"{sp}_goc", f"Gene-order conservation score (0-100) with {name}: how well the "
                                      f"neighbouring genes are conserved around the ortholog pair."))
        if sp in HAS_WGA:
            out.append((f"{sp}_wga", f"Whole-genome-alignment coverage (0-100) between the human and "
                                     f"{name} orthologs."))
    return out


SUMMARY_COLS = [
    ("n_species_with_ortholog", "Number of the 10 species with any ortholog call (0-10). Empty for unmapped rows."),
    ("n_species_one2one", "Number of the 10 species with a one-to-one ortholog (0-10)."),
    ("deepest_ortholog_species", "Most distant species in the ladder chimpanzee, macaque, mouse, rat, chicken, "
                                 "Xenopus, zebrafish, fly, worm, yeast that has an ortholog call; `none` if "
                                 "no species has one. A coarse gene-age proxy; a missing call is not proof "
                                 "of absence."),
]

TABLES = {
    "rbp_conservation_scores.tsv": (
        "Main table, one row per input RBP (2,335 rows).",
        [("Identity and transcript mapping", MAP_COLS + [SIM_COL]),
         ("CDS conservation (UCSC hg38 bigWigs)", track_cols() +
          [("cds_bases", "Number of CDS bases pooled for the scores above (sum of the CDS exon segments).")]),
         ("Annotations carried over from the input list", ANNOT_COLS),
         ("Ensembl Compara orthology, per species", ortho_cols()),
         ("Orthology summary", SUMMARY_COLS)]),
    "rbp_conservation_simple.tsv": (
        "Short version of the main table: one headline conservation score per RBP, with a reliability flag "
        "(built by `src/07_simple_table.py`).",
        [("Columns", [
            ("GeneSymbol", "Gene symbol from the input list."),
            ("UniProtID", "UniProt accession from the input list."),
            ("cds_phyloP_100way", "Headline score: mean phyloP over the CDS of the selected transcript, "
                                  "100-vertebrate alignment (= `phyloP100way_mean` in the main table). "
                                  "Roughly -0.6 to 6.3; higher = more conserved. Rough guides from the UCSC "
                                  "track notes: about 2 or more is typical of coding exons, about 5 or more is "
                                  "strongly constrained (18% of mapped RBPs reach 5). Empty for unmapped rows."),
            ("conservation_percentile", "Percentile rank (0-100) of `cds_phyloP_100way` among the mapped RBPs in "
                                        "this list; 100 = most conserved."),
            ("conservation_tier", "`high` = top quarter (percentile >= 75), `low` = bottom quarter (<= 25), "
                                  "`medium` otherwise. Relative to this list, not an absolute threshold. "
                                  "Empty for unmapped rows."),
            ("mouse_protein_identity_pct", "Percent of the human protein identical to the mouse ortholog "
                                           "(= `mmusculus_pid_query`). Empty if no mouse ortholog call."),
            ("n_species_with_ortholog", "Number of the 10 species (chimpanzee to yeast) with an ortholog call."),
            ("deepest_ortholog_species", "Most distant of those species with an ortholog call; `none` if no "
                                         "species has one."),
            ("score_reliable", "True if the RBP is mapped, the scored isoform contains >= 90% of the input "
                               "protein's residues, and >= 90% of CDS bases have alignment data."),
            ("reliability_note", "Why `score_reliable` is False; empty when True."),
        ])]),
    "rbp_transcript_map.tsv": (
        "RBP-to-transcript mapping and QC; the identity block of the main table.",
        [("Columns", MAP_COLS)]),
    "isoform_similarity_mismatches.tsv": (
        "Coverage of the input protein by the scored protein, for the 28 mapped RBPs without an exact-length match.",
        [("Columns", [
            ("GeneSymbol", "Gene symbol from the input list."),
            ("UniProtID", "UniProt accession from the input list."),
            ("input_aa", "Length of the input protein sequence (aa)."),
            ("scored_aa", "Length of the protein encoded by the scored transcript (aa)."),
            SIM_COL])]),
    "rbp_cds_segments.tsv": (
        "CDS exon segments of each selected transcript (Ensembl BioMart, GRCh38).",
        [("Columns", [
            ("ensembl_transcript_id", "Selected transcript."),
            ("exon_id", "Ensembl exon ID."),
            ("exon_rank", "Exon number within the transcript, 5' to 3' on the transcript."),
            ("chrom", "GRCh38 chromosome, `chr`-prefixed."),
            ("strand", "Transcript strand: 1 = forward, -1 = reverse."),
            ("cds_start_1based", "First coding base of this exon, 1-based, inclusive, genomic coordinate."),
            ("cds_end_1based", "Last coding base of this exon, 1-based, inclusive, genomic coordinate."),
            ("cds_len", "Coding bases in this exon segment (`end - start + 1`)."),
        ])]),
    "cds_conservation_by_exon.tsv": (
        "Conservation per CDS exon segment, for exon-level questions.",
        [("Columns", [
            ("ensembl_transcript_id", "Selected transcript."),
            ("exon_id", "Ensembl exon ID."),
            ("exon_rank", "Exon number within the transcript."),
            ("chrom", "GRCh38 chromosome, `chr`-prefixed."),
            ("start0", "Segment start, 0-based (UCSC convention)."),
            ("end", "Segment end, exclusive (`end - start0` = number of bases)."),
            ("n_bases", "Coding bases in the segment."),
        ] + [c for tr, label in TRACKS.items() for c in (
            (f"{tr}_mean", f"Mean {label} over the segment's covered bases."),
            (f"{tr}_covered", f"Number of segment bases with {label} data."))])]),
    "cds_conservation_by_transcript.tsv": (
        "Transcript-level conservation, the score block of the main table keyed by transcript ID.",
        [("Columns", [("ensembl_transcript_id", "Selected transcript.")] + track_cols() +
          [("cds_bases", "Number of CDS bases pooled.")])]),
}


def main():
    rows, md = [], ["# Column dictionary", "",
                    "What every column means in each table in `results/tables/`. Generated by "
                    "`src/06_column_dictionary.py`, which fails if a table has a column without an entry "
                    "here, so this file stays in step with the data. Machine-readable copy: "
                    "`results/tables/column_dictionary.tsv`.", "",
                    "Conventions: percentages are 0-100, fractions are 0-1, empty cells mean not applicable "
                    "or no data. Genome build GRCh38; Ensembl release 116; UCSC hg38 alignments.", ""]
    for fname, (blurb, sections) in TABLES.items():
        actual = pd.read_csv(T / fname, sep="\t", nrows=0).columns.tolist()
        documented = [c for _, cols in sections for c, _ in cols]
        missing = [c for c in actual if c not in documented]
        extra = [c for c in documented if c not in actual]
        assert not missing and not extra and len(documented) == len(set(documented)), \
            f"{fname}: undocumented={missing} not-in-table={extra}"
        md += [f"## `{fname}`", "", blurb + f" {len(actual)} columns.", ""]
        for title, cols in sections:
            if title != "Columns":
                md += [f"### {title}", ""]
            md += ["| Column | Meaning |", "|---|---|"]
            for c, d in cols:
                md.append(f"| `{c}` | {d.replace('|', '/')} |")
                rows.append({"table": fname, "section": title, "column": c, "description": d})
            md.append("")
    Path("COLUMNS.md").write_text("\n".join(md))
    pd.DataFrame(rows).to_csv(T / "column_dictionary.tsv", sep="\t", index=False)
    print(f"documented {len(rows)} columns across {len(TABLES)} tables")


if __name__ == "__main__":
    main()
