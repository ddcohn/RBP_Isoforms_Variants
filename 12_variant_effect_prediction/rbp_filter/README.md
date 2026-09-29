# Filtering the variant-effect tables down to RNA-binding proteins

The SpliceAI/DeepLoc/protGPS pipeline in `../splicing/` and
`../localization_condensate/` was run genome-wide -- on every GRCh38
ClinVar variant and every usable CMC variant, regardless of gene. That's
correct for the pipeline itself (nothing about SpliceAI, DeepLoc, or
protGPS is RBP-specific), but this project's actual scope is RNA-binding
proteins specifically, so the six merged tables need to be filtered down
before use, not read as-is.

## RBP gene list -- current version

`rbp_gene_symbols.txt` (**2,089 genes**) is the active filter list. It
comes entirely from UniProt itself (the `GN=` gene-name field on each
fetched protein record) -- **not** from either lab's isoform table.
Built from two starting ID sets:

1. **RBP2GO's own UniProt accessions** (1,977 IDs with `RBP_status ==
   RBP` in `High_confidence_human_RBPs_rbp2go.txt`, downloaded from
   [RBP2GO](https://rbp2gov2.dkfz.de/): References & Data -> Protein
   Data -> Homo sapiens -> RBP Dataset). An aggregation of many
   RNA-interactome-capture (CLIP-seq/mass-spec) studies, so it also
   captures RBPs with no annotated domain at all.
2. **A classic-RNA-binding-domain list** (534 NCBI gene IDs with
   `Has_RBD == 1` in *your own* `table_260823_with_rna.csv`, resolved to
   a `uniprot_accession` from that same table -- not the labmate's).

Both ID sets were fetched directly from UniProt's REST API
(`fetch_rbp_uniprot_sequences_v2.py`; all 2,095 accessions fetched
successfully, 0 failures), and the final gene list is simply the set of
distinct `GN=` values UniProt itself reports for those records
(`rebuild_clean_rbp_list.py`) -- 6 records had no `GN=` at all and were
dropped. Because the gene list is derived directly from records that
were already successfully fetched, **every one of these 2,089 genes has
a wild-type sequence by construction** (`rbp_wt_sequences.fasta`, not
checked into git -- large, and easily regenerated from
`rbp_uniprot_ids_to_fetch.txt`).

### A real bug, found and fixed

An earlier version of this list (2,215 genes, still visible in
`rbp_gene_symbols_PRECORRUPTED_backup.txt`) resolved RBP2GO's UniProt
IDs to gene symbols via the *labmate's* isoform table's `gene_symbol`
column instead of fetching from UniProt directly. That column turned
out to have a real data-quality problem: for 649 of the 1,968 RBP2GO
genes, the `gene_symbol` field itself contained a UniProt
evidence-code annotation glued onto the gene name, e.g.
`"AARS1 {ECO:0000303|PubMed:38653238, ECO:0000312|HGNC:HGNC:20}"`
instead of a clean `"AARS1"`. Since every downstream table/gene-symbol
column (ClinVar's `GeneSymbol`, CMC's `GENE_NAME`) uses clean names,
these 649 corrupted entries silently matched nothing -- every variant
for those 649 real genes was wrongly excluded from every RBP-filtered
table and notebook.

Fixing this (switching to UniProt's own `GN=` field, no lab table
involved) **net-gained 542 genes** (649 corrupted entries recovered,
minus 19 clean names that didn't reappear under this method --
mostly borderline/obscure calls, plus a few that look like they
shouldn't have been flagged as RBPs to begin with, e.g. `SOX10`,
`POU4F3`, `HMX3` are transcription factors, not RNA-binding proteins).
1,547 genes were confirmed correct as clean before the fix; the
corrupted 649 are the ones the fix actually recovers.
`rbp_gene_symbols_PRECORRUPTED_backup.txt`, `rbp_gene_symbols_rbp2go.txt`,
and `rbp_gene_symbols_domain_based.txt` are kept for provenance/history
but are **superseded** -- use `rbp_gene_symbols.txt`.

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

## Row counts, before -> after (current, corrected list)

| Table | Genome-wide | RBP-only (2,089 genes) |
|---|---|---|
| ClinVar DeepLoc deltas | 2,657,205 | 380,253 |
| ClinVar protGPS deltas | 2,587,382 | 363,750 |
| ClinVar SpliceAI scores | 4,137,146 | 574,335 |
| ClinVar metadata lookup | 4,494,430 | 645,025 |
| CMC DeepLoc deltas | 4,283,459 | 513,456 |
| CMC protGPS deltas | 4,190,344 | 494,651 |
| CMC SpliceAI scores | 5,070,804 | 596,266 |

(Two earlier passes exist in git history with smaller numbers: one on
the 533-gene domain-only list, one on the corrupted 2,215-gene list.
Both are superseded by the numbers above.)

## Wild-type sequences and saturation mutagenesis

`rbp_wt_sequences.fasta` / `rbp_wt_list.tsv` (2,089 sequences, 1,500,804
total residues -- the TSV has `GeneSymbol`, `UniProtID`, `Length`,
`Sequence` columns for direct use) is the basis for a full point-mutation
scan (every possible single-residue substitution at every position, run
through DeepLoc and protGPS) -- 19 x 1,500,804 = **28,515,276**
missense-only mutant sequences if generated for all 2,089 proteins. Not
yet built or run.

**DeepLoc and protGPS have been run on the 2,089 wild-type sequences
themselves** (the baseline the point-mutation deltas will be measured
against) -- `rbp_wt_deeploc_results.csv` and `rbp_wt_protgps_results.tsv`,
both verified by content: 2,089 rows each, exactly matching the input
count. protGPS additionally reports 8 sequences over its 5,000-residue
limit (`NA` placeholder rows) -- all genuinely giant proteins (AHNAK,
DST, EPPK1, KMT2D, MACF1, MDN1, SYNE1, SYNE2), not a bug. These are run
on the fresh UniProt sequences specifically, not reused from the
ClinVar/CMC pipeline's `WT_*` delta columns, since those come from a
different transcript source (RefSeq/Ensembl, not UniProt) and aren't
guaranteed to be the identical sequence or numbering.

## Notebooks

RBP-only versions of all six exploration notebooks are in
`../../notebooks/` (`*_rbp_exploration.ipynb`), reading the `*_rbp.tsv`
tables above. Re-executed against the corrected gene list.
