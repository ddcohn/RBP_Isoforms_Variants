# Evolutionary conservation of the RBP list

Conservation scores for the 2,335 RBPs in
`../12_variant_effect_prediction/rbp_filter/rbp_wt_list_with_go.tsv`
(`GeneSymbol`, `UniProtID`, ...), computed fresh from public sources
(consistent with the lab convention of recalculating rather than reusing
precomputed columns). Pulled 2026-10-09.

Two kinds of score, both gene-level, one row per RBP:

1. **Nucleotide conservation of the coding sequence** -- mean phyloP and
   phastCons over every CDS base of the selected transcript, from UCSC hg38
   multiple-species alignments (100 vertebrates; 30-way = 27 primates + 3
   outgroups).
2. **Protein/ortholog conservation** -- Ensembl Compara orthology (type,
   % identity, gene-order conservation, alignment coverage) against ten
   species from chimpanzee to yeast, plus counts and the most distant
   species with an ortholog call.

## Quick start: one conservation score per RBP

`results/tables/rbp_conservation_simple.tsv` has one row per RBP with a headline score
(`cds_phyloP_100way`, mean CDS phyloP over 100 vertebrates), its percentile and tier within
this list, mouse protein identity, ortholog depth, and a `score_reliable` flag. Filter on
`score_reliable == True` (2,321 of 2,335 RBPs). Everything below is the detailed table behind it.

```python
import pandas as pd
s = pd.read_csv("results/tables/rbp_conservation_simple.tsv", sep="\t")
s[s.score_reliable].sort_values("cds_phyloP_100way", ascending=False).head()
```

## Per-base and per-residue scores

`results/per_base/` (git-ignored, about 70 MB; rebuilt by steps 8 and 9 in roughly 10 minutes each):

- `rbp_cds_per_base_scores.parquet` -- phyloP and phastCons (100-way and 30-way) for every coding
  base of every selected transcript, 4,954,734 rows, with genomic position, exon, CDS position, codon
  number, codon position (1-3) and amino acid.
- `rbp_per_residue_scores.tsv.gz` -- one row per residue (1,649,252): codon mean and positions-1-and-2
  mean for each of the four tracks.

Example: `python src/10_phylop_example_figure.py HNRNPA1` draws per-base phyloP along one RBP's coding
sequence, coloured by codon position (`results/figures/HNRNPA1_phylop_per_base.png`). Every base is shown;
nothing is clipped by the axis.

Checks run by `src/09_per_base_scores.py`: the row count equals the total CDS length; each transcript's
per-base mean matches the step 3 mean to better than 1e-6; and codon position 3 has the lowest mean on both
strands (forward 0.90, reverse 0.78, against 4.9-5.4 at positions 1 and 2), which confirms the reading frame
and strand orientation. All 2,326 translations start with Met and end with a three-base stop codon.
Amino acids come from the Ensembl translation. For 6 RBPs the scored isoform has the input length but
different residues (HNRNPAB differs at 125 of 332 positions; LSM14A 6; ZNF598 4; MRPL45, MYEF2 and POLR2A 1
each), so their `aa` column does not match the input list's sequence.

New to the tables? `notebooks/conservation_walkthrough.ipynb` shows which file answers which question, with a
runnable example each (ranking RBPs, comparing groups, one gene per base, one gene per residue, a positions-1-and-2
score, joining your own annotation). It is saved with its outputs, so it displays without re-running.

## Main output

`results/tables/rbp_conservation_scores.tsv` -- 2,335 rows (9 unmapped rows
kept with `match_type == none` and empty scores), 93 columns.

| Column group | Columns |
|---|---|
| Identity / mapping | `GeneSymbol`, `UniProtID`, `ensembl_gene_id`, `ensembl_transcript_id`, `match_type` (uniprot / symbol / none), `selection_rule`, `symbol_agrees`, `length_matches_input`, `expanded_search`, `input_residues_in_scored_protein` |
| CDS conservation, per track `phyloP100way`, `phastCons100way`, `phyloP30way`, `phastCons30way` | `<track>_mean`, `_median`, `_frac_conserved` (phyloP >= 2, phastCons >= 0.9), `_cov` (fraction of CDS bases with alignment data) |
| Orthology, per species | `<sp>_orthology_type`, `<sp>_pid_query` (% of the human protein identical), `<sp>_pid_target`, `<sp>_n_orthologs`, `<sp>_goc`, `<sp>_wga` (the last two only where Ensembl provides them) |
| Summary | `n_species_with_ortholog` (of 10), `n_species_one2one`, `deepest_ortholog_species` |
| Carried over from the input list | `Ribosomal_protein`, `Spliceosomal_protein`, `Other_RNP_machine`, `HGNC_Gene_Family` |

**Every column in every table is described in [`COLUMNS.md`](COLUMNS.md)** (also as
`results/tables/column_dictionary.tsv`).

Other tables in `results/tables/`: `rbp_transcript_map.tsv` (mapping + QC),
`rbp_cds_segments.tsv` (CDS exon coordinates used), `cds_conservation_by_exon.tsv`
(per CDS exon segment, so exon-level questions need no re-query),
`cds_conservation_by_transcript.tsv`.
Figure: `results/figures/rbp_conservation_summary.{png,pdf}`.

## Method

- **Transcript choice** (`src/02_select_transcripts.py`): candidates are
  protein-coding transcripts annotated with the row's UniProt accession
  (fallback: gene symbol). Preference order: a transcript whose CDS encodes a
  protein of exactly the input `Length` (so the scored CDS is the same protein
  as the input sequence), then MANE Select, Ensembl canonical, longest CDS. If no
  candidate has the input length, the search is widened to every protein-coding
  transcript of the same gene plus the Ensembl transcripts cross-referenced from the
  UniProt entry, and those are used only if one matches the length exactly
  (`expanded_search == True`; 5 genes: ADAT3, ERCC6L2, PABPC1L, SMARCA1, TRIP12).
- **Length mismatches** (`src/02b_isoform_similarity.py`): for the 28 mapped RBPs
  where no Ensembl 116 transcript encodes the input length, the scored protein is
  compared with the input sequence. `input_residues_in_scored_protein` is the
  fraction of input residues found in matching blocks (difflib; a coverage measure,
  not an alignment identity). It is 1.0 for every exact-length match.
- **CDS conservation** (`src/03_cds_conservation.py`): BigWig values read over
  HTTP range requests (no multi-GB download), pooled over all CDS bases of the
  transcript so exons are length-weighted. Bases without alignment data are
  excluded, not zero-filled; `_cov` reports how many were excluded.
- **Orthology** (`src/01_fetch_biomart.py`, `src/04_merge_tables.py`): one
  Compara row per gene and species; one-to-one orthologs preferred, otherwise
  the highest query % identity.

## QC (this run)

| Check | Result |
|---|---|
| RBPs mapped to an Ensembl transcript | 2,326 / 2,335 (2,315 by UniProt accession, 11 by gene symbol) |
| Unmapped | EEF1A1P5, GPATCH4, HNRNPUL2-BSCL2, LOC127814297, MKRN4P, NXF5, RBMY1C, SNRPGP15, ZRSR2P1 |
| Selected transcript encodes a protein of the input length | 2,298 / 2,326 (98.8%); for the other 28 no Ensembl 116 transcript has the input length (isoform or annotation differences vs the UniProt sequence) |
| Mismatched RBPs that still contain >= 90% of the input residues | 23 / 28; below 90%: UNKL (0.38), PEG10 (0.46), DST (0.71), MACF1 (0.72), DDX31 (0.88) |
| Ensembl symbol equals input `GeneSymbol` | 2,320 / 2,326 (6 renamed genes, e.g. CNOT12 -> TNKS1BP1; UniProt accession is authoritative) |
| Summed CDS segments equal BioMart CDS length | 2,326 / 2,326 |
| CDS bases with phyloP100way data | min 99.7% |
| CDS bases with phyloP30way data | min 85.6% (one gene below 95%) |
| Spot check vs pyBigWig exact mean over the same bases | identical (2 exons) |
| Spot check vs the UCSC connector (4 exons) | within 0.09 phyloP (1-2%); the connector's own summarisation differs slightly |

Distribution across the 2,326 mapped RBPs (25th / median / 75th percentile):
phyloP100way 2.85 / 3.83 / 4.72; phastCons100way 0.76 / 0.87 / 0.93;
phyloP30way 0.64 / 0.76 / 0.86; phastCons30way 0.78 / 0.91 / 0.98.

## Caveats

- **dN/dS is not available.** Ensembl 116 BioMart exposes no `*_homolog_dn` /
  `*_homolog_ds` attributes for human, so dN/dS columns could not be pulled.
- **The 100-way and 30-way scores are on different scales.** The primate
  alignment is shallow, so gene-level phyloP30way means span roughly -0.4 to 1.1
  versus -0.6 to 6.3 for the vertebrate alignment. Compare within a track, not across.
- **Multi-copy and retrogene loci score low for alignment reasons**, e.g. the
  RBMY1* Y-chromosome family and RPS4Y1/2, EIF2S3B (high phyloP100way, near-zero
  phyloP30way). Treat these as unreliable.
- **No ortholog call is not proof of absence.** Five mapped genes have no
  Compara ortholog in any of the ten species (ATP5MF-PTCD1, FIP1L1, H4C1,
  NBPF10, RPL41); FIP1L1 and RPL41 are worth checking by hand.
- `deepest_ortholog_species` is a coarse gene-age proxy over ten species, not a
  phylostratigraphy.
- gnomAD gene constraint was not pulled: the bulk download host is blocked in
  this environment.

## Reproduce

```bash
conda env create -f environment.yml && conda activate rbp-conservation
cd 13_evolutionary_conservation
python src/01_fetch_biomart.py      # ~10 min, Ensembl BioMart release in data/raw/biomart_release.txt
python src/02_select_transcripts.py
python src/02b_isoform_similarity.py
python src/03_cds_conservation.py   # ~5 min, reads UCSC bigWigs remotely
python src/04_merge_tables.py
python src/05_figures.py
python src/07_simple_table.py
python src/08_fetch_peptides.py     # ~5 min, Ensembl BioMart
python src/09_per_base_scores.py    # ~10 min, per-base and per-residue tables
python src/10_phylop_example_figure.py HNRNPA1   # example per-base figure for one RBP
python src/06_column_dictionary.py   # rebuilds COLUMNS.md; fails if a column is undocumented
```

`data/raw/` (BioMart dumps, ~58 MB) is git-ignored and re-fetched by step 1.
Ensembl BioMart answers GET requests only and sometimes returns an HTML
"Service unavailable" page with status 200; `src/biomart.py` detects this and
retries.
