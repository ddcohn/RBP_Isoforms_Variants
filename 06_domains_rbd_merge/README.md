# Domain annotations and classical RBD flags

Rather than compute domain/disorder-region annotations from scratch, this
step pulls an already-computed, richer version from elsewhere in the lab's
shared project space and merges it onto the clean base table — after
verifying it's actually trustworthy first.

## Where this data came from, and why we checked it first

A labmate's file, `/u/project/kappel/RBP/Misc/Classical_RBDs_Annotated.csv`
(2.9GB), already has both:
- `Domains`, `Domains_count`, `Domains_avg_size`, `Domains_total_size`,
  `Domains_range`, `Domains_discrete_seq` — per-domain-family annotation
  (source: InterPro), keyed by protein.
- `Has_RBD`, `RBD_name_ranges`, `RBD_all_ranges`, `RBD_names` — a
  classical-RNA-binding-domain flag and named ranges.

Before merging it in, we checked two things that have bitten this project
before:

1. **`check_rbd_domains.py`** — checks for the same kind of duplicate-row
   problem found in the original broken isoform table. Result: 25,589 rows
   collapse to 16,951 distinct proteins, but every duplicate row for a given
   protein carries byte-identical data (0 inconsistencies) — safe to dedupe
   by taking any one row per protein. Also surfaces the real coverage
   number: only 16,951 of this project's 20,431 proteins (83%) have an
   entry in this file at all.

2. **`check_domain_overlap.py`** — the more interesting check. `Domains_range`
   and `Domains_discrete_seq` are dicts keyed by domain-family name (e.g.
   `{'RRM': [(55,134),(135,218),(230,303)]}`), and the concern was whether
   overlapping domain calls (common with InterPro, since multiple domain
   databases can annotate the same region under different names) get
   merged or double-counted. **Confirmed they are not merged**: for one
   concrete case (`Q92499`), summing `Domains_total_size` across families
   gives 794 residues for a 740-residue protein, because a 178-residue
   B30.2/SPRY domain sits entirely nested inside a 427-residue Helicase
   ATP-binding domain call, and both are counted independently. Real
   (non-boundary-touching) overlaps affect 113 proteins within the same
   domain family and 222 proteins across different families — small
   fractions of the 8,181 proteins with any domain annotation, but real.
   **Practical implication**: fine to use the per-family data as-is; do not
   sum `Domains_total_size` across families for a "total domain coverage"
   number without first merging overlapping ranges.

3. **`merge_domains_rbd.py`** — the actual merge, same left-join +
   backup-first pattern as every other merge in this project. 16,950/20,431
   rows (83.0%) matched.

4. **`verify_domains_merge.py`** — standard post-merge integrity check.
   Confirms 20,431 rows in both files, all 27 pre-existing columns
   byte-identical, no corruption.

## A related side-investigation: `extract_idr_from_fraza_table.py`

Before finding the RBD file above, a labmate mentioned "what Faris had"
(the Hoffman2 username `fraza`) should be fine for IDR data specifically.
This script extracts the much richer IDR characterization
(`idr_method`, `IDR_FCR`, `IDR_NCPR`, `IDR_kappa`, `IDR_delta`, etc. — a
full per-IDR biophysical profile, not just count/size/range) from the
*original broken* isoform table
(`Isoform_Post_Merge_PSLab_OpenTargets_Updated_Interim_20260817.csv`),
after confirming the duplicate rows in that file are internally consistent
for these specific columns (same check pattern as above). Kept here as a
record of that investigation; superseded in practice by the
`Classical_RBDs_Annotated.csv` merge above once the richer domain+RBD file
was found, but the extracted IDR data itself is real and dedupe-verified if
needed later.

## Column reference

| Column | Meaning |
|---|---|
| `Domains` | list of domain-family names present |
| `Domains_count` / `_avg_size` / `_total_size` | dict keyed by family name (see overlap caveat above) |
| `Domains_range` / `_discrete_seq` | dict of family → list of `(start,end)` ranges / subsequences |
| `Has_RBD` | protein-level flag: does this protein have any classical RNA-binding domain |
| `RBD_name_ranges` / `RBD_all_ranges` / `RBD_names` | named classical-RBD ranges |
