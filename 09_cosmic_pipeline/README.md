# COSMIC integration

Adds somatic (cancer) mutation data to complement ClinVar's germline
classifications. **Important framing correction made partway through this
work**: this project's isoform table is not RBP-restricted — it's
proteome-scale (20,431 rows; insulin, serum albumin, and hemoglobin are
all in it). "RBP" in every directory/file name here refers to the lab's
research focus, not a filter on the table's contents. Every coverage
number below is against the full table, not an RBP subset.

## Why COSMIC, and what it is (and isn't)

ClinVar answers "is this inherited (germline) variant known to cause
disease" with a curated Pathogenic/Benign/VUS classification. COSMIC
answers a different question: "how often does *cancer* mutate this gene,
and where" — from raw somatic mutation observations in sequenced tumor
samples. **COSMIC's raw mutation files are not filtered for pathogenicity
or driver status** — confirmed by checking that a `Synonymous_Rows` count
(silent mutations, which by definition can't be damaging) is present and
non-trivial in every gene's data. Most mutations in a sequenced tumor are
incidental "passengers," not disease-driving.

COSMIC also splits mutation data across separate files by scope —
`Cosmic_MutantCensus` is Cancer-Gene-Census-restricted (~700-750 genes
only), while `Cosmic_GenomeScreensMutant` is genome-wide. Given the table
is proteome-scale, the genome-wide file is the one that actually matters
here.

## The join key: ENSG, not gene symbol

COSMIC's gene reference file (`Cosmic_Genes`) links `COSMIC_GENE_ID` ↔
`GENE_ACCESSION` (ENSG, versioned — strip the suffix) ↔ `ENTREZ_ID` ↔
`HGNC_ID` ↔ `GENE_SYMBOL`, all in one row. The mutation files themselves
only carry `GENE_SYMBOL` + `COSMIC_GENE_ID`, so `Cosmic_Genes` is the
bridge table. **Deliberately avoided gene-symbol matching** throughout —
this project has already been bitten by symbol-matching fragility
elsewhere (the lab's `isRBP` classification, ClinVar's shifting HGNC
symbols) — ENSG is an exact match with no ambiguity.

## Scripts

1. **`build_cosmic_gene_summary.py`** — streams the 62.77M-row
   `Cosmic_GenomeScreensMutant` file (gzipped, ~20 minutes at observed
   throughput), aggregating per `COSMIC_GENE_ID`: total mutation rows,
   distinct genomic mutations (`GENOMIC_MUTATION_ID` — the real dedup key,
   equivalent role to ClinVar's `VariationID`), distinct samples, and a
   breakdown by mutation type (missense/truncating/synonymous) and
   somatic-status confidence. Joins in the `Cosmic_Genes` reference fields.
   Output: `Cosmic_Gene_Summary_v104_GRCh38.csv`, 19,509 genes.
   Deduplicating by `GENOMIC_MUTATION_ID` shrinks 62.77M row-observations
   to 17.52M distinct mutations (72.1% reduction; average 3.58 samples per
   distinct mutation) — the legitimate ClinVar-style dedup, as opposed to
   the row-*duplication* problems found elsewhere in this project.

2. **`check_clinvar_ensembl_id_sample.py`** / **`check_clinvar_ensembl_id_full.py`**
   — before assuming ClinVar's own `ENSEMBL_ID` column (needed to bridge
   ClinVar ↔ COSMIC ↔ this table, all on ENSG) is reliable, checked its
   real population rate and internal consistency. A 10K-row sample
   suggested only 91.2% coverage; the full 114.9M-row scan (submitted via
   `submit_ensembl_check.sh`, ~72 minutes) showed the sample was
   unrepresentative — **real coverage is 97.79%**, with near-perfect
   internal consistency (a handful of genuine mismatches, e.g. `BRCA1`
   spuriously co-occurring with `SYNE2`/`APC`'s real ENSGs and one row
   with literal JSON-fragment text where a gene symbol should be — traced
   to rare CSV field-misalignment from unescaped characters in the
   `Germline_Class` column, not a systemic problem).

3. **`check_cmc_rbp_table_coverage.py`** — after finding the labmate had a
   local copy of COSMIC's **Cancer Mutation Census** (CMC) file
   (`cmc_export.tsv` / `CancerMutationCensus_AllData_v104_GRCh37.tsv.gz`,
   which — despite the name — turned out to cover ~19,996 genes
   genome-wide, not just Cancer-Gene-Census genes as initially assumed;
   confirmed by direct check before believing it), this checks real
   coverage and per-column population against the isoform table's
   transcript list. Result: 70.0% of CMC's 5.8M rows match one of the
   table's transcripts by exact `ACCESSION_NUMBER` (ENST); only 68.8% of
   the table's transcripts appear in CMC at all — likely because exact
   transcript-ID matching misses real matches where COSMIC and this
   project's `exact_match_rebuild.py` picked *different* (but same-gene)
   transcripts. **Recommendation for the actual merge, not yet done**:
   match on ENSG (via `ACCESSION_NUMBER` → Ensembl lookup), not exact
   transcript ID, to recover that gap.

   The CMC file's extra columns are genuinely useful but unevenly
   populated: `MIN_SIFT_SCORE` 100%, `GERP++_RS` 72.3%, `GNOMAD_EXOMES_AF`
   15.7%, `GNOMAD_GENOMES_AF` 13.2%, `CLINVAR_CLNSIG` only 1.8%,
   `DNDS_DISEASE_QVAL_SIG` 0.03% (deliberately significance-filtered, so
   expected to be rare). Low population for gnomAD/ClinVar isn't
   necessarily a gap — most somatic-only mutations were never
   independently submitted to ClinVar as germline variants, and gnomAD
   only reports variants it actually observed.

## Not yet built

- The ENSG-based (rather than transcript-ID-based) join between CMC and
  the isoform table.
- A mutation-level (not gene-level) table analogous to the ClinVar
  `variant_summary.txt` pull in `../10_clinvar_repull/` — scoped but not
  yet executed; requires the same `GENOMIC_MUTATION_ID` dedup logic
  described above.
