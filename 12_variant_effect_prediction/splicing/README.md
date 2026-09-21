# Splice-effect prediction (SpliceAI)

Runs Illumina's SpliceAI (v1.3.1) on ClinVar and COSMIC/CMC variants to
predict splicing consequences. Unlike localization/condensate, this needs
no protein-sequence construction — SpliceAI takes genomic position +
REF/ALT directly, and ClinVar/COSMIC already provide VCF-ready alleles.

## ClinVar (`build_clinvar_spliceai_chunks.py`, `submit_clinvar_spliceai_array.sh`)

Builds a VCF from all 4,494,430 GRCh38 ClinVar rows (`chr`-prefixed
contigs matching the reference FASTA; `MT` remapped to `chrM`; 3 variants
on unplaced (`Un`) scaffolds excluded — no matching contig in the primary
assembly). SpliceAI has no internal batching — one forward pass per
variant, ~150-250ms each depending on node contention — so a serial run
would take multiple days. Split into 200 chunks (~22,473 variants each)
and run as a parallel SGE array job on the general compute pool instead
of a single dedicated queue, to maximize how many chunks can run at once.

**First attempt at this failed almost completely** — see the top-level
`12_variant_effect_prediction/README.md` for the full diagnosis
(undersized `h_rt`/`h_data` based on a small-scale test, plus SpliceAI
not writing output incrementally so a killed task loses everything).
Fixed by rerunning with `h_rt=6:00:00`, `h_data=24G`, and a per-task log
path.

## COSMIC/CMC (`build_cmc_spliceai_chunks.py`, `submit_cmc_spliceai_array.sh`, `count_unique_cosmic_variants.py`)

COSMIC's full `GenomeScreensMutant` file is 62.7M rows, but almost all of
that is redundancy — SpliceAI's prediction depends only on genomic
position + alleles, not which tumor sample carried it, and the same
variant recurs across many samples.
`count_unique_cosmic_variants.py` confirmed this: only **15,645,235**
unique `(chrom,pos,ref,alt)` combinations across the full 62.7M rows.

Used the smaller, separately-provided Cancer Mutation Census (CMC) file
instead — 5,802,908 rows, **5,756,934 unique** `GENOMIC_MUTATION_ID`s
(<1% redundancy to begin with), and critically, it carries *both*
GRCh37 and GRCh38 genome-position columns natively, so no coordinate
conversion was needed to reuse the existing GRCh38 reference FASTA.
After dropping ~50K rows with a missing/unmapped GRCh38 position and
deduping by `(chrom,pos,ref,alt)`: **5,491,494 usable unique variants**,
split into 260 chunks (~21,122 variants each), same array-job approach
and same corrected resource limits as the ClinVar run from the start.

## `merge_spliceai_results.py`

Concatenates the per-chunk SpliceAI output VCFs and extracts the four
delta scores (`DS_AG`/`DS_AL`/`DS_DG`/`DS_DL` — acceptor/donor gain/loss
probability, 0-1) and four delta positions (`DP_*`, bp offset from the
variant) per variant, writing a standalone `VariationID -> scores` table
— **deliberately kept separate from the main ClinVar table**, not
merged into it, per lab preference to keep prediction outputs modular.
SpliceAI's own published guidance: scores >=0.2 suggest a possible
splicing effect, >=0.5 is high-confidence; most variants score ~0.00
across all four (no predicted effect).
