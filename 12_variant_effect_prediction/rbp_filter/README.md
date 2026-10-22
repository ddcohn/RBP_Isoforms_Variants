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

`add_rnp_machine_flags.py` adds `Ribosomal_protein`,
`Ribosomal_protein_terms`, `Spliceosomal_protein`,
`Spliceosomal_protein_terms`, `Other_RNP_machine`, and
`Other_RNP_machine_terms` columns, matched against each gene's
`GO_Cellular_Component` entries.

**Rewritten from an earlier substring-keyword version that had real
false positives.** The original matched on whether a keyword like
"ribosome" or "snrnp" appeared *anywhere* in the GO_Cellular_Component
text -- that caught unrelated terms that happen to share a word
fragment: 42 proteins flagged `Ribosomal_protein=Y` via "preribosome"
(an assembly *intermediate*, i.e. a ribosome biogenesis factor, not the
mature ribosome -- `BMS1`, `BOP1`, `NOP14`, etc.); 4 proteins flagged
`Spliceosomal_protein=Y` via "U7 snRNP" (histone mRNA 3'-end
processing) or "7SK snRNP" (RNA Pol II elongation control) -- real
RNPs, just unrelated to splicing (`HEXIM1`, `LARP7`, `LSM11`, `MEPCE`);
3 proteins flagged `Other_RNP_machine=Y` via "signal recognition
particle *receptor* complex," which has no RNA component at all and
isn't an RNP (`RP9`, `SRPRA`, `SRPRB`).

Fixed with exact GO ID allowlists (9 IDs for the mature ribosome, 26
for the spliceosome) instead of substring keywords -- every term was
pulled from this table's own real data and individually verified
against its actual GO definition, not guessed from memory. The U7/7SK
snRNP and preribosome hits are still surfaced in `Other_RNP_machine`
under their correct label (they're real RNPs, just the wrong specific
category) rather than silently dropped; the SRP receptor proteins were
dropped entirely since they genuinely aren't an RNP. The new
`*_terms` columns record the exact matched GO ID behind every `Y`, so
any individual flag can be audited directly against the source data.
`audits/rnp_machine_flags_audit.tsv` lists all 479 flagged proteins
with their matched term, for spot-checking without re-deriving anything.

Verified with an automated consistency check across all 2,335 proteins
(every `Y` has a matching term, every `N` has none -- 0 failures), plus
a cross-check against independent HGNC gene families (the
`RNP_machine_flags_cleanup_README.md`-equivalent summary): confirmed 0
false negatives for the ribosomal family, but found a real, separate
gap -- 15 proteins HGNC curates as core spliceosome-complex members
(`CCAR1`, `CHERP`, `FAM32A`, `FAM50A`, `HNRNPAB`, `NOSIP`, `PPIG`,
`PQBP1`, `PRCC`, `PUF60`, `RBM10`, `RBM25`, `SAP18`, `THRAP3`, `WBP11`)
that have no spliceosome-specific GO_Cellular_Component term at all --
a GO annotation-coverage gap, not a bug in this script. Not yet
recovered (would need the same GO_Biological_Process/HGNC_Gene_Family
pattern used for the RNA-class flags below).

Current counts (of 2,335 genes): **203 ribosomal, 194 spliceosomal, 97
other named RNP machine** -- down from the original (buggy) 245/198/100.
The "other" category covers named RNP machines beyond those two --
exosome (RNase complex), signal recognition particle, telomerase
holoenzyme complex, box C/D and H/ACA snoRNP, the spliceosome commitment
complex (verified redundant with real spliceosomal GO terms for all 5
proteins it touches, so harmless), editosome, vault RNP, nuclear pore
complex. Spot-checked true positives against known members (`RPS6`/
`RPL7` ribosomal, `SF3B1`/`SNRNP70` spliceosomal, `EXOSC1-10` exosome,
`SRP9/14/19/54/68/72` signal recognition particle, `TERT`/`TEP1`/`DKC1`
telomerase) and they still check out.

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

## RNA target type: mRNA / tRNA / rRNA / snRNA / snoRNA / piRNA / miRNA / siRNA / lncRNA binding

`add_trna_rrna_flags.py` (superseded, kept for history -- see below)
originally added just `tRNA_binding`/`rRNA_binding`, matching only the
single base GO term for each. `add_rna_class_flags.py` replaces both
with **18 columns: a Y/N flag plus a `_terms` column per class**, for
9 RNA classes. Each flag matches *every* real GO ID for that class
found in this table's own `GO_Molecular_Function` data (not assumed
from memory -- found by grepping the column for "RNA binding" and
manually sorting ~45 distinct real hits into classes); the matching
`_terms` column records the actual `GO_ID:term_name` entry/entries
that matched (semicolon-joined if more than one sub-term hits, e.g.
`COIL` shows both `U1 snRNA binding` and `U2 snRNA binding` under
`snRNA_binding_terms`), same pattern as `Other_RNP_machine_terms`
above. The script is idempotent -- it clears any columns it (or its
superseded predecessor) previously wrote before adding fresh ones, so
re-running never produces duplicate columns (a real bug caught and
fixed while first adding the `_terms` columns: a naive "drop the old
two column names" filter didn't account for the broader set of 9
class names a prior run had already written).

| Column | GO IDs matched | Count |
|---|---|---|
| `mRNA_binding` | GO:0003729, GO:0036002, GO:1990825 | 253 |
| `tRNA_binding` | GO:0000049, GO:1904678, GO:1990856 | 63 |
| `rRNA_binding` | GO:0019843, GO:0008097, GO:0070180, GO:0070181, GO:1990400 | 60 |
| `snRNA_binding` | GO:0017069, GO:0017070, GO:0030619-0030626, GO:0071209, GO:0097322 | 41 |
| `snoRNA_binding` | GO:0030515, GO:0034511, GO:0034512, GO:0034513 | 23 |
| `piRNA_binding` | GO:0034584 | 6 |
| `miRNA_binding` | GO:0035198, GO:0070878, GO:0070883 | 39 |
| `siRNA_binding` | GO:0035197 | 7 |
| `lncRNA_binding` | GO:0106222, GO:0002153 | 17 |

Broadening tRNA/rRNA to their full sub-term families changed their
counts from the original 62/39 to 63/60 -- rRNA picked up 21 more genes
that only had a ribosomal-subunit-specific term (`large`/`small
ribosomal subunit rRNA binding`) rather than the generic `rRNA binding`
term. `mRNA_binding` previously had no flag at all (documented as "the
implicit default for neither tRNA nor rRNA set") -- now explicit, since
`GO:0003729` (mRNA binding) exists in the data and was simply unused.

Deliberately excluded: structural/motif terms that don't identify an
RNA class (double-/single-stranded, G-quadruplex, poly(U)/poly(C),
regulatory region, generic "translation factor activity, RNA binding"),
and single-specific-ncRNA-gene terms already covered by
`Other_RNP_machine` (`ribonuclease P RNA binding`, `7S RNA binding` =
SRP RNA) or not cleanly assignable to one class (`telomerase RNA
binding`, `telomeric repeat-containing RNA binding` -- TERC is its own
thing, not folded into lncRNA here). Spot-checked: `piRNA_binding`
catches all 4 PIWIL genes (PIWIL1-4) plus 2 others with real GO
annotations; `lncRNA_binding` catches AGO2, DICER1, FMR1, MECP2,
TARBP2, PRKRA -- all genuinely lncRNA-pathway-associated.

Still true regardless of class count: the mRNA **sub-region** half
(5'UTR/CDS/intron/3'UTR) is the separate ENCORI-derived table,
`encori_binding_site_summary_rbp_filtered.tsv`, which has no RNA-species
column of its own -- the two tables answer different halves of the
original "RNA target type" question and aren't merged into each other.

**Sanity-checking the original two flags caught a real, pre-existing bug in this file:**
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

### Recovering coverage via GO_Biological_Process and HGNC_Gene_Family

The table above (`GO_Molecular_Function` only) covers 450 of 2,335
proteins (19.3%) with a specific RNA class, well under the 1,286 (55.1%)
that carry the general `GO:0003723` "RNA binding" MF term -- most
RNA-binders have the generic term but no class-specific one.
`extract_bp_candidates.py` pulled every real `GO_Biological_Process`
term and `HGNC_Gene_Family` value across the whole table that contains
each RNA-class keyword (e.g. "mRNA", "tRNA"), for manual curation
against real data rather than guessing which BP terms imply binding.

`add_bp_family_rna_class.py` applies the curated result: only for the
969 proteins that have the general RNA-binding term but no specific
class (the "gap"), check BP terms that imply processing/binding that
RNA type (e.g. `GO:0006397` mRNA processing, `GO:0006418` tRNA
aminoacylation) and HGNC families (e.g. "Aminoacyl tRNA synthetases,
Class I/II" -> tRNA). Explicitly excluded "regulation of X
transcription by RNA polymerase Y" BP terms (regulates the DNA-template
transcription step, not binding of the mature RNA) and one real false
positive caught during curation: HGNC family "MicroRNA protein coding
host genes" (159 proteins) means the gene's own intron happens to host
a miRNA gene -- it says nothing about whether the protein itself binds
miRNA.

Result: 366 of the 969 gap proteins (37.8%) recovered a specific class,
raising total specific-class coverage from 19.3% to **34.9%** (816 of
2,335). Biggest gains: mRNA 253->481, rRNA 60->181 (nearly 3x), tRNA
63->95. 603 proteins still have only the generic term with no specific
class recoverable this way -- closing that further would need either
the full GOA annotation file (this table's GO columns are likely a
filtered subset) or ENCORI's experimental CLIP data (which RNA species
a protein was actually pulled down with, a stronger signal than any GO
term). `build_rna_binding_plot.py` / `rna_binding_breakdown.py`
produce the coverage bar chart. **Not yet merged back into
`rbp_wt_list_with_go.tsv` itself** -- this recovery currently only
exists as the standalone analysis/plot, not as new columns in the main
table.

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

### Also joined onto a labmate's much larger merged table

A labmate (asharma) independently built a far more extensive merged
table for the same 2,332 RBPs --
`/u/project/kappel/asharma/RBP/merged/merged_260903B_ptm_pred_slim_RBPs_domains_spurs_rsa_iso_hpa.csv`
(137 columns: real + predicted PTMs, SLiMs, SPURS structural
ddG/pLDDT, RSA, a richer IDR biophysical profile than this repo's
own, HPA subcellular localization, condensate predictions, plus
everything this repo's own isoform table has, including
`HGNC_Gene_Family`). Checked column-by-column before assuming
anything was missing: the only genuine gaps were this repo's 12
RNP-machine/RNA-class *base* columns (`Ribosomal_protein`,
`Spliceosomal_protein`, `Other_RNP_machine`, and the 9 RNA-class
binding flags) -- none of their `_terms` companions existed either.

`merge_into_asharma_table.py` joins those 22 columns (12 base + 10
`_terms`) onto asharma's table by `uniprot_accession`, writing to
**the user's own directory, not asharma's** (the source table is
read-only input, never modified in place):
`/u/project/kappel/ddcohn/RBP/merged_260903B_ptm_pred_slim_RBPs_domains_spurs_rsa_iso_hpa_with_rnp_rna_class.csv`
(159 columns, all 2,332 rows matched, 0 gaps). Not checked into git
(293 MB).

`add_protgps_deeploc_wt.py` then adds this project's own WT-baseline
protGPS (12 condensate scores) and DeepLoc (17 localization/signal
scores) predictions onto the same file, prefixed `protgps_`/`deeploc_`,
joined by UniProt accession -- 188 columns total, still 2,332 rows, 0
gaps. TANGO has no WT-baseline run to add (only ever run on the 31.5M
saturation mutants) -- can be added the same way later if that baseline
gets run.

`add_encori_clip.py` then adds ENCORI's CLIP-seq binding data too,
joined by `gene_symbol` -- the one exception to "one row per protein"
in this whole merge chain. Per explicit instruction not to average
across cell lines, the 79 of 195 ENCORI-matched genes with more than
one cell-line row get their entire row duplicated once per ENCORI
match rather than collapsed (e.g. `CSTF2T` -> 4 rows, one each for
HEK293/HeLa/HepG2/K562). Table grows from 2,332 to 2,459 rows as a
result. Anything joining on `gene_symbol`/`uniprot_accession`
downstream should expect duplicates for multi-cell-line genes and
decide deliberately how to handle them, rather than assuming one row
per protein still holds after this step.

`extract_encori_motifs.py` / `add_encori_motifs.py` add ENCORI's
"Motifs" sheet (top 2 enriched motifs per RBP x CLIP-method x cell-
line experiment, 196 of our RBPs matched, 375 rows) the same no-
averaging way -- but as a genuine **cross-product**, not a compound
join against the per-cell-line rows `add_encori_clip.py` already
added: a motif's own `cellTissue` field is itself a pipe-separated
list of multiple cell lines (e.g.
`BE(2)-C|HEK293|HeLa|HepG2|K562`), not a single value that could
align 1:1 with one binding-site row. So a gene with N existing rows
and M motif rows becomes N*M rows -- motifs are gene-level findings,
not tied to one specific binding-site cell line, so this isn't
incoherent, just real further growth. Table grows from 2,459 to 2,758
rows. By this point `gene_symbol`/`uniprot_accession` duplication is
the norm for well-studied genes, not the exception -- treat this
merged table as one-row-per-(protein, experiment-combination), not
one-row-per-protein, for anything downstream.

`add_nls_nes_sps.py` adds 32 more columns from a later version of
asharma's table
(`merged_260903B_ptm_pred_slim_RBPs_domains_spurs_rsa_iso_hpa_nls_nes_sps.csv`):
classical NLS predictions (`cNLS_*`, mono- and bi-partite), NES
predictions (`NES_*`), and a second independent method for both
(`SPS_NLS_*`/`SPS_NES_*`). Checked first that this source's other 137
columns were byte-identical to what `merge_into_asharma_table.py`
already pulled in -- confirmed, only these 32 are new. Unlike the two
ENCORI scripts above, this is a **per-protein broadcast join, not a
row multiplier**: NLS/NES/SPS data is one set of values per
`uniprot_accession`, so it's just repeated across whichever duplicate
rows a protein already has from the ENCORI merges, not a new source
of row growth. Table stays at 2,758 rows, now 253 columns.

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
