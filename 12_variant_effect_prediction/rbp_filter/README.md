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
`rbp_uniprot_ids_to_fetch.txt`). `build_rbp_wt_list.py` converts that
FASTA into `rbp_wt_list.tsv` (`GeneSymbol`, `UniProtID`, `Length`,
`Sequence`, one row per gene, sorted by symbol) -- this is the file
every later `add_*.py` script in this directory actually reads from,
not the FASTA directly.

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

## Row counts, before -> after (current list, 2,335 genes)

| Table | Genome-wide | RBP-only (2,335 genes) |
|---|---|---|
| ClinVar DeepLoc deltas | 2,657,205 | 424,735 |
| ClinVar protGPS deltas | 2,587,382 | 408,232 |
| ClinVar SpliceAI scores | 4,137,146 | 648,255 |
| ClinVar metadata lookup | 4,494,430 | 726,288 |
| CMC DeepLoc deltas | 4,283,459 | 569,249 |
| CMC protGPS deltas | 4,190,344 | 550,444 |
| CMC SpliceAI scores | 5,070,804 | 662,812 |

(Three earlier passes exist in git history with smaller numbers: the
533-gene domain-only list, the corrupted 2,215-gene list, and the
corrected-but-not-yet-combined 2,089-gene list. All superseded by the
numbers above.)

## Combined with a user-supplied gene list

`protein_ids_user_supplied.txt` (1,392 unique gene symbols) was
cross-referenced against the 2,089-gene list above: 1,097 already
overlapped, 295 didn't. Of those 295, `combine_lists.py` resolved 246 as
genuinely new genes (fetched fresh from UniProt, same method as
above), 45 as already covered under the current UniProt-preferred gene
name (e.g. an older symbol like `AARS` for what's now `AARS1`), and 4 as
unresolvable -- 2 are Ensembl/GenBank clone identifiers rather than
real gene symbols (`AC013461.1`, `RP1-37E16.12`), 1 is a
gene-gene readthrough/fusion name whose two component genes (`RBM14`,
`RBM4`) were already individually present, and 1 (`PRPF4B`) turned out
to already be covered too, just under UniProt's own primary name for
that record (`PRP4K`) rather than the current HGNC symbol.
`build_combined_list.py` merges the 246 new genes in.

`protein_ids_user_supplied.txt` is from a mass spec dataset of unknown
experiment type/QC history (not confirmed whether contaminants or a
no-bait/no-RNA control were already filtered out). Several of the 246
newly-added genes are not RNA-binding proteins by known function (e.g.
`CAT`, `APOB`, `HLA-A`, `CRYAB`, `AFP`, `KRT18`, proteasome subunits,
TCA-cycle enzymes) -- checked the 246 against the cRAP database (the
standard proteomics contaminant reference, ~94-124 proteins from
[Zenodo](https://zenodo.org/records/15115102)): 6 are confirmed
contaminants by gene symbol or UniProt accession -- `CAT`, `CYCS`,
`GSTP1`, `KRT18`, `NQO1`, `UBE2I`. Decision: keep all 246 in the active
list anyway, including those 6 -- the cRAP check was informational, not
treated as grounds for exclusion, given the dataset's QC history is
unknown either way. Worth revisiting if the mass spec dataset's
provenance/QC becomes known later.

**Current active gene list: 2,335 genes** (`rbp_gene_symbols.txt`,
`rbp_wt_sequences.fasta` / `rbp_wt_list.tsv`, 1,658,067 total residues).
All six `*_rbp.tsv` tables and all six RBP-only notebooks have been
regenerated against this combined list (see the updated row-count table
below and note in the notebooks).

## Wild-type sequences and saturation mutagenesis

The wild-type sequence set above is the basis for a full point-mutation
scan (every possible single-residue substitution at every position, run
through DeepLoc and protGPS) -- 19 x 1,658,067 = **31,503,273**
missense-only mutant sequences if generated for all 2,335 proteins. Not
yet built or run.

**DeepLoc and protGPS have been run on all 2,335 wild-type sequences
themselves** (the baseline the point-mutation deltas will be measured
against) -- `rbp_wt_deeploc_results.csv` and `rbp_wt_protgps_results.tsv`
(the original 2,089-gene run plus an incremental run on the 246 new
genes, merged), both verified by content: 2,335 rows each, exactly
matching the gene count. protGPS additionally reports 8 sequences over
its 5,000-residue limit (`NA` placeholder rows) -- all genuinely giant
proteins (AHNAK, DST, EPPK1, KMT2D, MACF1, MDN1, SYNE1, SYNE2), not a
bug. These are run on the fresh UniProt sequences specifically, not
reused from the ClinVar/CMC pipeline's `WT_*` delta columns, since those
come from a different transcript source (RefSeq/Ensembl, not UniProt)
and aren't guaranteed to be the identical sequence or numbering.

## Notebooks

RBP-only versions of all six exploration notebooks are in
`../../notebooks/` (`*_rbp_exploration.ipynb`), reading the `*_rbp.tsv`
tables above. Re-executed against the corrected gene list.

## Gene Ontology annotation

`rbp_wt_list_with_go.tsv` (`add_go_terms.py`) adds `GO_Cellular_Component`,
`GO_Biological_Process`, `GO_Molecular_Function` columns to `rbp_wt_list.tsv`,
pulled from *your own* `table_260823_with_rna.csv` and joined directly by
UniProt accession (both tables use UniProt IDs natively, so no gene-symbol
translation needed). 2,332 of 2,335 genes matched (2,321 with actual
non-empty GO data); the 3 unmatched are genuinely obscure identifiers not
expected to be in a standard gene-centric table (`ATP5MF-PTCD1` and
`HNRNPUL2-BSCL2` are gene-readthrough/fusion names, `LOC127814297` is an
unnamed provisional NCBI locus).

`table_260823_with_rna.csv` itself is not checked into this repo (it's an
external, pre-existing file from earlier pipeline work) and reproducing
this step doesn't require it -- the join's output, `rbp_wt_list_with_go.tsv`,
already has the GO columns baked in and is the thing committed here.

(Checked the labmate's PSLab isoform table too, as an alternative source --
it has cleaner `role_in_transcription`/`role_in_translation`/
`role_in_mrna_stability` flags, but doesn't track splicing as a category
at all: a known splicing factor, U2AF1, has all three flags set to 0. Used
the GO-term source instead since it actually captures splicing-related
function.)

## RNP machine membership

`add_rnp_machine_flags.py` adds `Ribosomal_protein`, `Spliceosomal_protein`,
`Other_RNP_machine`, and `Other_RNP_machine_terms` columns, derived by
keyword-matching each gene's `GO_Cellular_Component` string (e.g. a hit
on "cytosolic ribosome" / "ribosomal subunit" -> `Ribosomal_protein=Y`;
"spliceosomal complex" / "U1 snRNP" etc. -> `Spliceosomal_protein=Y`).
The "other" category covers named RNP machines beyond those two --
exosome (RNase complex), signal recognition particle, telomerase
holoenzyme complex, box C/D and H/ACA snoRNP, the spliceosome commitment
complex, editosome, vault RNP, nuclear pore complex.

245 ribosomal, 198 spliceosomal, 51 other named RNP machine (of 2,335
genes). This is a straightforward substring match over real GO terms,
not independent curation -- it inherits whatever GO's own annotation
scope and quality is for each gene, and a keyword can be genuinely
ambiguous. Caught one during review: an earlier "signalosome" keyword
matched `APC`/`GRB2`/`LCP2`, which are unrelated signal-transduction
proteins (the COP9 signalosome, not an RNA-related complex) -- removed
before finalizing. Spot-checked the rest against known members (`RPS6`/
`RPL7` ribosomal, `SF3B1`/`SNRNP70` spliceosomal, `EXOSC1-10` exosome,
`SRP9/14/19/54/68/72` signal recognition particle, `TERT`/`TEP1`/`DKC1`
telomerase) and they check out.

## RBP family (hnRNP, SR protein, DEAD-box helicase, etc.)

`fetch_hgnc_families.py` pulls the official HGNC (HUGO Gene Nomenclature
Committee) curated gene-family group for each of the 2,335 genes from
their REST API (`https://rest.genenames.org/fetch/symbol/<SYMBOL>`,
`gene_group` field) -- this is the authoritative source for named RBP
families (e.g. HGNC's own "Heterogeneous nuclear ribonucleoproteins",
"DEAD-box helicases", "Serine and arginine rich splicing factors" groups
are literally what "hnRNP"/"DEAD-box helicase"/"SR protein" mean).
`add_hgnc_family.py` joins the result in as a new `HGNC_Gene_Family`
column in `rbp_wt_list_with_go.tsv` (semicolon-separated if a gene
belongs to more than one HGNC group).

2,032 of 2,335 genes have at least one HGNC family assigned; 295 matched
an HGNC record with no family group; 8 had no HGNC match at all.
Spot-checked well-known examples directly from the output file:
`HNRNPA1` -> Heterogeneous nuclear ribonucleoproteins, `DDX3X` -> DEAD-box
helicases, `SRSF1` -> Serine and arginine rich splicing factors, `AGO2`
-> Argonaute RISC component family, `ELAVL1` -> ELAV like RNA binding
protein family -- all correct.

## RNA target type: tRNA / rRNA binding

`add_trna_rrna_flags.py` adds `tRNA_binding` and `rRNA_binding` Y/N columns,
derived by checking `GO_Molecular_Function` for two specific GO terms:
`GO:0000049` (tRNA binding) and `GO:0019843` (rRNA binding) -- confirmed
against this table's own data (`grep`-checked the term names attached to
those IDs) before trusting them. 62 genes flagged tRNA-binding, 39
rRNA-binding, no overlap between the two sets. This is the mRNA-vs-
rRNA-vs-tRNA half of the "RNA target type" question asked earlier in the
project (the mRNA sub-region half -- 5'UTR/CDS/intron/3'UTR -- is the
separate ENCORI-derived table, `encori_binding_site_summary_rbp_filtered.tsv`,
which has no RNA-species column of its own). Between the two, mRNA-binding
itself has no explicit flag -- it's the implicit default for any RBP with
neither tRNA_binding nor rRNA_binding set, not a separately confirmed GO
match.

**Sanity-checking these flags caught a real, pre-existing bug in this file:**
spot-checking `rRNA_binding` (then the last column) via `awk -F'\t'
'$14=="Y"'` returned nothing, despite the script reporting 39 matches and
a raw value count confirming 39 `Y`s were actually in the column. The file
had Windows-style CRLF line endings on every line (`file` confirms, and
`git show` on the much earlier `ccfa034` commit shows the same -- this
predates today's change, likely present since that commit or earlier),
so the true last column's value on every line was `Y\r`/`N\r`, not
`Y`/`N` -- identical on screen, but failing exact string match in `awk`,
`grep -x`, etc. Python's `csv` module handles CRLF correctly regardless
of position, so this didn't silently corrupt any of this project's actual
joins/merges (all done via `csv.DictReader`/`DictWriter`) -- it only broke
naive shell-level exact-match checks on whichever column happened to be
physically last. Fixed by stripping all `\r` bytes (`tr -d '\r'`) across
the whole file. Worth remembering: **when a shell-based spot-check on a
TSV/CSV's last column doesn't match what the data clearly shows, check
for CRLF before assuming the logic is wrong.**

### Also joined onto the full isoform table, RBP-filtered

The isoform table built in `01_isoform_table_analysis/` through
`09_cosmic_pipeline/` (`table_260823_with_rna.csv`, 20,431 proteins) is
**not** RBP-specific -- this project's 2,335-RBP scope only exists in
`rbp_wt_list_with_go.tsv` above, derived *from* the isoform table via a
`Has_RBD` filter, not the reverse. Isoform-level analyses that need the
isoform table's other ~35 columns (PPI partners, domains, GO/tissue/
OpenTargets, COSMIC, etc.) alongside RBP family classification needed
their own row-filtered copy rather than cramming everything into the
variant-effect gene list.

`build_rbp_isoform_table.py` filters `table_260823_with_rna.csv` down
to the 2,335 RBPs (by `uniprot_accession`) and joins in
`HGNC_Gene_Family` from the fetch above, writing
`/u/project/kappel/ddcohn/RBP/isoform_table_rbp_only_with_hgnc_family.csv`
(2,332 of 2,335 matched -- the same 3 obscure/readthrough IDs that don't
match anywhere in this pipeline: `LOC127814297`, `ATP5MF-PTCD1`,
`HNRNPUL2-BSCL2`). Not checked into git (71 MB, and this repo doesn't
check in the isoform table itself either -- see its own README).

## RNA target region (mRNA sub-region: 5' UTR / CDS / intron / 3' UTR)

Not derived from this pipeline's own data -- sourced from the user's
own ENCORI reference workbook (`ENCORI_human_RBP_reference.xlsx`,
built from Zhou et al. 2026 *Nature Methods*, CLIP-seq-derived), a
separate large ENCORI binding-site processing effort outside this
repo. Its "Binding Site Summary" sheet gives, per RBP per cell
line/tissue, the % of that RBP's binding regions falling in each
functional gene region (3'-UTR, 5'-UTR, CDS, intron, coding/non-coding
exon, splice sites, start/stop codon, plus several repeat-element
categories). **This is region-within-gene-model, not RNA species** --
it doesn't have an mRNA-vs-rRNA-vs-tRNA classification column, only
answers the "where in an mRNA" half of the original ask.

Coverage is partial by nature of CLIP-seq data: 195 of this project's
2,335 RBPs are in the sheet at all (ENCORI only has experimental data
for RBPs that have actually been CLIP-seq-profiled). `extract_encori_
binding_summary.py` filters the sheet to just those RBPs and writes
`encori_binding_site_summary_rbp_filtered.tsv` -- **row-for-row as-is,
no averaging across cell lines**: 79 of the 195 matched RBPs have more
than one cell-line/tissue row (e.g. CSTF2T has 4), and collapsing them
would discard real per-experiment binding data. So this is a one-to-
many table (RBP -> multiple rows), not a one-row-per-RBP table like
the isoform table or `rbp_wt_list_with_go.tsv` -- join on `RBP` and
expect duplicates for any RBP profiled in more than one cell line.
