# RBP Isoform Table: RNA/CDS, STRING PPI, and OpenTargets Enrichment

Scripts used to diagnose a broken RBP isoform table, rebuild it from a clean
base table, and enrich it with RNA/CDS sequences, protein-protein interaction
data (STRING), and disease association data (OpenTargets).

All scripts were run against files on the Hoffman2 HPC cluster
(`/u/project/kappel/...`) via SSH, except where noted. Most were tested on a
small `--limit`/`--dry-run` sample before being run at full scale, and use
JSON checkpoint files (paths hardcoded near the top of each script) so a
crashed or killed run can resume without redoing completed work.

## Background

The original table (`/u/project/kappel/RBP/Isoform_Table/Isoform_Post_Merge_PSLab_OpenTargets_Updated_Interim_20260817.csv`,
316 columns, 68,645 rows) had serious structural problems: 3,872 proteins with
more than one row flagged as the "dominant" isoform, 29,502 dominant-flagged
rows for only 20,483 distinct proteins, and 534 proteins with none at all. The
`UNIQUE` key column was 98% empty and had duplicates where populated. See
`01_isoform_table_analysis/`.

Rather than patch this, we started from a clean base table provided by a
labmate (`table_260823.csv`, 20,431 rows, one row per `uniprot_accession`, no
duplicates) and built everything else on top of it.

## Pipeline order

1. **`01_isoform_table_analysis/`** — diagnostic scripts run against the
   *old* broken table only. Not part of the rebuild; kept for reference.
2. **`02_rna_cds_pipeline/`** — pulls the full mRNA (`rna_sequence`) and
   coding sequence (`cds_sequence`) for each protein.
3. **`03_sanity_checks/`** — validates the RNA/CDS columns by translating
   each CDS with the standard genetic code and comparing it to the protein
   sequence already in the table.
4. **`04_string_ppi_pipeline/`** — pulls STRING interaction partners
   (`PPI_UniProt_Partners`, `PPI_UniProt_Partners_in_Dataframe`), computed
   three independent ways for cross-validation.

Each subdirectory has its own README with per-script details.

## Phase 2: further enrichment, and auditing the wider lab's related data

5. **`05_go_tissue_opentargets_enrichment/`** — adds GO terms, bulk tissue
   expression, and OpenTargets disease associations.
6. **`06_domains_rbd_merge/`** — merges in domain annotations and classical
   RNA-binding-domain flags from a labmate's file, after verifying it's
   trustworthy first (duplicate-row and domain-overlap checks).
7. **`07_clinvar_variant_audit/`** — diagnostic pass over a labmate's
   existing ClinVar variant-classification table; found a ~10x count
   inflation bug and traced its likely cause.
8. **`08_missing_rna_cds_investigation/`** — categorizes *why* the 876
   RNA/CDS-less proteins from stage 2 lack a sequence, then verifies that
   categorization against Ensembl's authoritative data (and finds the
   first-pass heuristic was directionally right but individually
   unreliable).
9. **`09_cosmic_pipeline/`** — adds COSMIC somatic mutation data, with the
   ENSG-based (not gene-symbol) join pattern established here reused in
   later stages.
10. **`10_clinvar_repull/`** — re-pulls ClinVar directly from NCBI after
    finding real problems in the derived file audited in stage 7, rather
    than trying to patch it.
11. **`11_column_inventory/`** — generic per-column stats tooling used
    throughout.
12. **`12_variant_effect_prediction/`** — scores ClinVar/COSMIC variants
    as mutant-vs-wild-type effects across splicing (SpliceAI) and
    subcellular localization / condensate propensity (DeepLoc, protGPS).
    Kept as standalone tables (see stage 10's decision to keep ClinVar
    and COSMIC separate), not merged into the isoform table. In progress
    — splicing and localization/condensate have working pipelines;
    PTM gain/loss and protein stability are not started.

**See [`FINDINGS.md`](FINDINGS.md) for the cross-cutting discoveries** —
including that this table is proteome-wide, not RBP-restricted, and that
at least three separate, inconsistent isoform-table lineages exist across
the lab.

## Why RNA/CDS coverage will never reach 100%

Even though the base table has a correct `ensembl_gene` assigned per protein,
and `exact_match_rebuild.py` exhaustively checks *every* transcript of that
gene (not just a guessed "canonical" one) for an exact translation match,
three genuine gaps remain:

1. **Some rows have no `ensembl_gene` at all.** 278 proteins in the base
   table have neither `ensembl_gene` nor `ensembl_protein` populated — there
   is no gene ID to start the search from. This is a gap in the source data,
   not something fixable by querying Ensembl differently.
2. **A gene can have an ENSG but no transcript matching the exact protein.**
   Once every transcript of a gene has been checked and translated, and none
   of them match, that is a verified negative (Ensembl's current annotation
   for that locus simply doesn't include a transcript producing this exact
   protein sequence) — not a sign the pipeline didn't try hard enough.
3. **Some proteins have no single fixed mRNA sequence to find at all.**
   Immunoglobulin/T-cell-receptor chains are the clearest case: the mature
   protein is assembled by somatic DNA recombination independently in each
   immune cell, so there is no one "correct" reference mRNA sequence sitting
   in the genome — see `02_rna_cds_pipeline/categorize_no_map.py` for the
   full breakdown of this and the other biological categories (pseudogenes,
   lncRNA-hosted micropeptides, mitochondrial micropeptides) that behave the
   same way.

## Key findings worth knowing before reusing this code

- **STRING silently caps results at 10 partners/protein** unless you pass an
  explicit `limit` override — the default response looks complete but isn't.
- **STRING silently renames/normalizes some input IDs** before returning
  results, so query IDs must be resolved via `get_string_ids` first and
  matched back using the *resolved* ID, not the original input.
- **A gene ID (ENSG) is not a protein ID.** One gene can have dozens of
  transcripts producing different proteins; naively picking "the canonical
  transcript" for a gene matches the wrong protein most of the time (measured
  60-91% wrong on this dataset). The reliable method translates every
  transcript of a gene and keeps only the one that exactly matches the stored
  protein sequence (`02_rna_cds_pipeline/exact_match_rebuild.py`).
- **Selenoproteins will "fail" a naive CDS-to-protein check** even when
  correct: the codon TGA is recoded to selenocysteine (U) in ~25 human genes,
  which a standard genetic-code translator reads as a stop codon. Check
  against the actual list of human selenoprotein genes before treating these
  as errors.
- **UniProt's protein sequences are not always byte-identical to a naive
  translation even when correct** — e.g. this table's base protein sequences
  have selenocysteine (`U`) replaced with cysteine (`C`), a deliberate
  standardization for tool compatibility, not corruption.
- **The Hoffman2 login node kills background processes** that exceed a ~1GB
  memory / 1hr CPU ulimit, silently and without a traceback. Long-running
  jobs must be submitted via `qsub` (UGE/SGE), not run directly with
  `nohup ... &` on the login node. See `02_rna_cds_pipeline/submit_exact_match.sh`
  for a working job template (`h_data`, `h_rt` resource requests).
- **STRING's live REST API serves an older release (v12.0)** than the newest
  downloadable STRING file (v12.5) — check `version.string-db.org` before
  assuming the API is current. A version-pinned endpoint
  (`version-12-5.string-db.org`) exists and works despite its own
  `/version` endpoint misreporting itself as 12.0.

## Requirements

Plain-stdlib Python 3 (tested against the system Python 3.6.8 on Hoffman2 —
no `requests`, `pandas`, etc. available; all HTTP calls use `urllib`).
