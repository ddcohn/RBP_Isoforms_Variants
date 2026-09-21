# Building mutant sequences from ClinVar's HGVS notation

Localization/condensate predictors (DeepLoc, protGPS) need actual protein
sequences, but ClinVar only gives HGVS notation describing the change
(e.g. `NM_014630.3(ZNF592):c.3136G>A (p.Gly1046Arg)`). This builds
wild-type and mutant sequence pairs directly from that notation, without
relying on any pre-existing lab table for the wild-type sequence (not
independently verifiable, so not trusted here — see repo-level
convention).

## `parse_clinvar_protein_changes.py`

Regex-parses the `(p....)` portion of ClinVar's `Name` field for every
GRCh38 row and classifies it. Order matters: deletion/duplication/
insertion patterns (`p.Glu315del`, `p.Thr123dup`, `p.Ala42_Ala43insAspAla`)
must be checked *before* the generic single-substitution pattern
(`p.Gly1046Arg`) — a naive substitution regex will happily match
`"Glu315del"` as if `"del"` were a 3-letter amino acid code. This bug was
caught by comparing category counts before/after: fixing it recovered
~20,000 variants that had been silently dropped into `unparsed`.

Category counts across all 4,494,430 GRCh38 variants:

| Category | Count | Directly editable from WT sequence? |
|---|---|---|
| missense | 2,516,726 | yes |
| no protein notation (noncoding/intronic) | 917,960 | n/a — no protein effect |
| synonymous | 764,214 | n/a — no protein effect |
| frameshift | 151,543 | **no** — needs CDS translation |
| nonsense | 103,437 | yes (truncate) |
| del | 21,993 | yes |
| ins | 6,452 | yes |
| delins | 5,936 | yes |
| dup | 3,020 | yes |
| unparsed (rare notations, e.g. repeat-expansion `p.138PG[3]`) | 1,966 | no |
| stoploss ("nonstop extension") | 1,183 | **no** — needs 3'UTR sequence |

`frameshift` and `stoploss` both look editable at first glance but
aren't: a frameshift's abbreviated notation (`p.Leu473fs`) only says
*where* the frame breaks, not what the new downstream sequence reads as
— that requires translating the shifted reading frame from the actual
CDS. `stoploss` means the stop codon itself is mutated into an amino
acid, so translation runs on into the 3' UTR — sequence that was never
part of the annotated protein to begin with. Both are left for a
possible later CDS-based pass; everything else (2,657,205 variants) is
handled here directly.

COSMIC/CMC has the equivalent breakdown for free, as an existing
categorical column (`Mutation Description AA`) rather than free-text
notation needing a regex parser — see `../splicing/count_unique_cosmic_variants.py`
and the lab notebook for those counts. Not yet wired into this
sequence-construction pipeline.

## `fetch_wt_sequences.py`

For the ~19,400 unique RefSeq transcript accessions (`NM_...`) referenced
by the directly-editable variants, fetches each transcript's GenBank
record from NCBI (`efetch`, batched ~200 accessions/request to stay
polite without an API key) and pulls the exact wild-type protein
sequence straight from the CDS feature's `/translation=` qualifier —
already present in the same record, no separate protein-database lookup
or translation step needed. 19,405/19,405 fetched successfully.

## `build_mutant_sequences.py`

Applies the parsed edit directly to the fetched WT sequence:
missense/nonsense are a one-residue swap or truncation; del/ins/delins/dup
edit the relevant range. Validates that the WT residue at the claimed
position actually matches the fetched sequence before accepting an edit
— catches transcript/numbering mismatches rather than silently producing
a wrong sequence. Mismatch rate: 359/2,657,205 (0.014%), and inspecting
those individually confirmed they're not data problems: mostly
delins/del ranges that extend through the stop codon itself (genuinely
unrepresentable without the CDS, correctly rejected) and HGVS's own
`Xaa` ("unspecified residue") notation, plus one mitochondrial-genome
accession (`NC_012920.1`, multi-gene record) that isn't a single-transcript
accession.

Output: 2,657,205 mutant sequences + 19,405 unique WT sequences, as FASTA,
plus a `VariationID -> mutant sequence ID` mapping table for joining
predictions back to variants afterward.
