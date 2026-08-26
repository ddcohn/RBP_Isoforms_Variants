# Sanity checks

The core validation method used throughout: translate `cds_sequence` with the
standard genetic code and compare it to the `sequence` (protein) column
already in the table. An exact match is strong evidence the right transcript
was selected; a mismatch needs explaining.

- **`sanity_check_seqs.py`** — first full-table pass: exact-match rate,
  structural CDS checks (starts `ATG`, length divisible by 3, ends in a stop
  codon), and whether `cds_sequence` appears as a substring of
  `rna_sequence`. Found 96.7% exact match on the first RNA pull, with
  mismatches clustering into an "X + protein" pattern (translated sequence
  is the real protein with one extra leading residue) traced to
  immunoglobulin/TCR partial gene-segment transcripts that start mid-codon.

- **`deep_check.py`** — digs into the worst single mismatch found
  (`A0A087WUL8`, translated 57 aa vs a real 3,843 aa protein) and traces it
  to a root cause: the protein's `ensembl_gene` field listed **two different
  genes**, and the transcript-selection logic grabbed a transcript from the
  wrong one. Also checks whether the benign "X-prefix" pattern explains the
  rest of the structural-check failures (it doesn't, mostly).

- **`test_multigene_hypothesis.py`** — quantifies how much of the
  `deep_check.py` finding generalizes: are the remaining ~628 unexplained
  mismatches enriched for proteins mapped to multiple genes? Yes (13.7% vs
  8.3% baseline), but it only explains 86 of 628 — most of the pile needed
  more digging.

- **`investigate_542.py`** / **`investigate_542_v2.py`** — categorize the
  remaining 542 single-gene mismatches by length pattern (translated
  shorter/longer/same-length) and real UniProt gene names/descriptions
  (v1 failed silently because the working table has no gene-name columns;
  v2 fetches names fresh from UniProt). Found that ~49% of "translated
  shorter" cases stop exactly at a `TGA` codon — the signature of
  selenocysteine recoding — and confirmed against the full list of ~25 known
  human selenoprotein genes (`GPX1-4/6`, `TXNRD1-3`, `DIO1-3`, `SELENO*`,
  `SEPHS2`): all 24 selenoproteins in the mismatch set are explained this
  way. Also directly confirmed (by comparing against live UniProt) that the
  table's protein sequences replace selenocysteine (`U`) with cysteine (`C`)
  — a deliberate standardization, not corruption, but it means even a
  selenocysteine-aware translator wouldn't match the stored sequence at that
  position.

- **`sanity_check_gene_fallback.py`** — same exact-match methodology applied
  specifically to the 156 rows recovered by `gene_level_fallback.py`. Found
  only 3.8% exact match / 91% unexplained mismatch, and (via a
  same-vs-different-prefix breakdown) that 64% diverge from the very first
  residue — i.e. the wrong transcript entirely, not a real biological
  variant. This result is what motivated `exact_match_rebuild.py`.

- **`full_table_sanity_check.py`** — the definitive whole-table pass,
  breaking down exact-match rate by `transcript_selection_method` with known
  explainable categories (X-prefix, selenoprotein, multi-gene) factored out.
  This is the script to rerun after any pipeline change to get a current,
  authoritative picture — it's what revealed `fallback_no_canonical_tag`
  (401 rows) had the same low-reliability problem as `gene_level_fallback`.
