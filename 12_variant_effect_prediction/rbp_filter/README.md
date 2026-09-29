# Filtering the variant-effect tables down to RNA-binding proteins

The SpliceAI/DeepLoc/protGPS pipeline in `../splicing/` and
`../localization_condensate/` was run genome-wide -- on every GRCh38
ClinVar variant and every usable CMC variant, regardless of gene. That's
correct for the pipeline itself (nothing about SpliceAI, DeepLoc, or
protGPS is RBP-specific), but this project's actual scope is RNA-binding
proteins specifically, so the six merged tables need to be filtered down
before use, not read as-is.

## RBP gene list

The RBP definition used is `Has_RBD == 1` in
`/u/project/kappel/ddcohn/RBP/table_260823_with_rna.csv` (this project's
own gene-level isoform table, 20,431 genes) -- 534 NCBI gene IDs, the
narrowest/most conservative definition (a confirmed RNA-binding domain),
not a broader RBP census list.

That table has no gene-symbol column (only `ncbi_geneid`,
`ensembl_gene`, `uniprot_accession`), so `build_rbp_gene_list.py` maps
those 534 NCBI gene IDs to gene symbols via the raw ClinVar file's own
`GeneID`/`GeneSymbol` columns (which cover ~20,256 genes). 533 of 534
resolved; one gene ID (`100913187`) has no match in ClinVar's gene list
and is not reflected in the filtered tables below -- a known, negligible
(1-gene) gap, not silently dropped without a note.

Output: `rbp_gene_symbols.txt`, the 533 resolved gene symbols, checked
into this repo since it's small and is the actual filter definition used
-- anyone rerunning `filter_rbp_tables.py` reproduces the same filter
from it directly, without needing to recompute the ID resolution step.

## Filtering the six tables

`filter_rbp_tables.py` filters each of the six merged tables in
`/u/project/kappel/ddcohn/protein_variant_effects/` down to rows whose
gene is in `rbp_gene_symbols.txt`, writing a `*_rbp.tsv` sibling next to
each original (not checked into git -- covered by the same data
`.gitignore` rule as the unfiltered tables). DeepLoc/protGPS deltas
already carry a gene column directly; ClinVar's SpliceAI table and CMC's
SpliceAI table don't (SpliceAI needs no gene/protein info to run), so
those are filtered via a join: ClinVar through
`clinvar_metadata_lookup.tsv` (`VariationID` -> `GeneSymbol`, built from
every GRCh38 row), and CMC through a fresh `GENOMIC_MUTATION_ID` ->
`GENE_NAME` map read directly from the raw CMC file (`cmc_protein_changes.tsv`
only covers the protein-editable variant categories, not the
synonymous/frameshift/noncoding variants SpliceAI also scored, so it
isn't a complete enough join key on its own).

## Row counts, before -> after

| Table | Genome-wide | RBP-only (533 genes) |
|---|---|---|
| ClinVar DeepLoc deltas | 2,657,205 | 108,379 |
| ClinVar protGPS deltas | 2,587,382 | 108,379 |
| ClinVar SpliceAI scores | 4,137,146 | 152,069 |
| ClinVar metadata lookup | 4,494,430 | 174,837 |
| CMC DeepLoc deltas | 4,283,459 | 150,661 |
| CMC protGPS deltas | 4,190,344 | 150,661 |
| CMC SpliceAI scores | 5,070,804 | 176,519 |

All 533 RBP genes are represented by at least one wild-type sequence
somewhere in the pipeline (531/533 in ClinVar, 524/533 in CMC, union =
533/533) -- none were dropped for lack of a buildable sequence.

## What's not done yet

The notebooks in `../../notebooks/` still read the unfiltered,
genome-wide tables. RBP-only versions of the exploration notebooks
haven't been built yet.
