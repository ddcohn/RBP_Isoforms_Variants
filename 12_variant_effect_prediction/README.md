# Variant effect prediction

Five new prediction categories were scoped for scoring ClinVar/COSMIC
variants as **mutant-vs-wild-type effects**, not just wild-type protein
annotation: splicing, subcellular localization, condensate propensity, PTM
gain/loss, and protein stability. Per-lab convention (see
`../07_clinvar_variant_audit/`), everything here is recalculated
independently rather than reusing labmates' precomputed columns.

Only two of the five categories have working code so far — **splicing**
(`splicing/`) and **localization/condensate propensity**
(`localization_condensate/`, via DeepLoc + protGPS). PTM gain/loss
(MusiteDeep) and protein stability (ThermoMPNN) are not started.

## Why mutant-vs-WT needs its own sequence-construction step

ClinVar and COSMIC don't contain protein sequences — only HGVS notation
describing the change (e.g. `p.Gly1046Arg`, `c.3136G>A`). Localization and
condensate-propensity predictors need actual sequences in and out. Rather
than trust a pre-existing lab table for wild-type sequences (provenance
not verifiable), `mutant_sequence_construction/` fetches WT sequences
fresh from NCBI per-transcript and edits them directly using the HGVS
protein notation already in ClinVar's `Name` field. See that folder's
README for the full method and the category breakdown (missense/nonsense/
frameshift/etc.), which determines what's directly usable this way versus
what would need full CDS-level translation (deferred).

Splicing doesn't have this problem — SpliceAI takes genomic
position + REF/ALT directly (which ClinVar/COSMIC do have), no protein
sequence needed.

## Compute environment

All large jobs run on Hoffman2 via `qsub` (UGE/SGE), never directly on the
login node — same ~1GB/1hr ulimit issue noted in the top-level README, but
it bit harder here: naive `torch`/`tensorflow` imports with GPU libraries
can exceed the limit even for a quick import check, with no traceback.
Job scripts need `#!/bin/bash -l` (login shell) or the `module` command
silently fails ("module: command not found").

**GPU submission on Hoffman2 requires resource flags beyond `cuda=1`.**
The Kappel lab's own dedicated GPU node (`kappel_gpu.q`, 2x L40S) forces
`highp` and `require_gpu=1` as job-side requirements — a job without both
sits in the queue (`qw`) indefinitely with no error. The general shared
GPU queues (`gpu_a100.q`, `gpu_l40s.q`, `gpu_h100.q`, `gpu_RTX2080Ti.q`,
plus others discovered only by broad request: `gpu_smp.q`, `gpu_rh7.q`,
`gpu_v100.q`) force `require_gpu=1` but *not* `highp` — requesting
`highp` restricts a job to `kappel_gpu.q` only, while omitting it opens
the job to the whole shared pool. SGE does correctly set
`CUDA_VISIBLE_DEVICES` to the granted device (verified directly), so
jobs land on an actually-idle GPU rather than fighting other users for a
busy one.

## Reproduce from scratch

The whole splicing + localization/condensate pipeline, for both ClinVar
and CMC, is now a Snakemake workflow (`workflow/Snakefile`,
`workflow/config.yaml`) instead of manually running each stage's script
and `qsub` job in order. `workflow/config.yaml` is the single source of
truth for every path/chunk-count/conda-env name the scripts used to
hardcode independently.

**This now also covers the satmut + SLiM/ELM + ClinVar/COSMIC benchmarking
pipeline** (`saturation_mutagenesis/README.md`), the RNP-machine-flags fix,
and the RNA-class GO-coverage recovery (`rbp_filter/README.md`) -- all
wrapped into the same `Snakefile`/`config.yaml`, not a second parallel
system. One command reproduces everything in this directory, given the
external inputs below are in place.

**One-time setup on Hoffman2:**
1. Create the `spliceai`, `deeploc2`, `protgps`, and `snakemake` conda
   envs from the exported specs in `../environments/` -- see that
   folder's README for exact commands and what each env is for.
2. Place the raw input files at the paths in `workflow/config.yaml`:
   `ClinVar_variant_summary_complete.csv`,
   `CancerMutationCensus_AllData_v104_GRCh37.tsv.gz`, the reference
   genome FASTA, the TANGO binary (third-party, see
   `saturation_mutagenesis/README.md` "TANGO setup"), and a labmate's
   `rbp_satmut_all.csv.gz` (satmut source; see that README's "Data
   source" section -- this repo cannot generate that file itself).
3. **Deploy the code to Hoffman2**: `./deploy_to_hoffman2.sh` from this
   directory, locally. The Snakefile's paths resolve relative to
   wherever it lives on disk, not to this git repo directly -- Hoffman2
   needs its own synced copy. Re-run this after editing any script.

**Run it:**
```
module load miniforge/23.11.0 && source activate snakemake
cd /u/project/kappel/ddcohn/RBP_pipeline_scripts/workflow   # the deployed copy, not the git checkout
snakemake --profile profiles/hoffman2 --cores 1 all
```
This submits every stage as its own `qsub` job via the
`profiles/hoffman2/` cluster profile (capped at 60 concurrent jobs via
`jobs:` in that profile -- keeps the account well under Hoffman2's
concurrent-job limit even though the old array-job approach queued
hundreds of tasks under one job ID), watches each job via
`qacct`/`qstat` instead of a human polling `qstat`, and automatically
retries a failed chunk with escalated memory/runtime (the `retries:` +
attempt-scaled `resources:` in the Snakefile encode this session's
actual OOM-kill-then-bigger-`h_data` lesson directly, instead of a
one-off hand-written retry script each time). Every prediction rule
ends with a count check against its own input and fails the job (not
just logs a warning) on a mismatch -- a truncated/killed run is a real
Snakemake failure, generalizing this session's `input_variants ==
output_variants` / exit-code-propagation fixes into the standard
pattern for every rule instead of something bolted on after the fact.
`snakemake -n --cores 1 all` dry-runs the whole DAG without submitting
anything (config.yaml is loaded automatically by the Snakefile itself
-- don't pass `--configfile` on the command line too, it greedily
consumes the next argument and will misparse the target as a filename).

Target `all` produces the 6 merged delta/score tables (SpliceAI/DeepLoc/
protGPS x ClinVar/CMC), `clinvar_metadata_lookup.tsv`,
`satmut_all_predictions.tsv` and its downstream delta/ClinVar/COSMIC
analyses, and the fixed RNP-machine flags -- everything in
`/u/project/kappel/ddcohn/protein_variant_effects/`.

**A real, known Snakemake gotcha, hit and fixed while extending this**:
after adding the satmut rules, a dry run against the real already-
completed outputs reported almost everything as stale and wanting to
rebuild (678 chunk-level prediction jobs!) -- not because anything was
actually wrong, but because some *inputs* (e.g. the big merged table,
which kept getting new columns added throughout this project after SLiM
disruption was first computed from it) have a newer file-modification
time than their *outputs*, and Snakemake's staleness check is mtime-
based, not content-based. Fixed the one time this happens with
`snakemake --touch --cores 1 all`, which marks every real, already-
verified-correct output as up to date without re-running anything --
confirmed immediately after with a dry run reporting "Nothing to be
done." If you ever see a suspiciously large rebuild plan against outputs
you know are already correct, suspect this before assuming the DAG
logic is wrong.

**Verified** (not just written): a dry run after the satmut extension
resolves the exact real job count (693 total: the original 7 clinvar/cmc
targets + 226x3 satmut chunk predictions + every downstream analysis
rule) and, after `--touch`, reports the already-completed real pipeline
as fully up to date; a live round-trip through the actual SGE profile
(force-rerunning `rnp_machine_flags`, a fast, safe-to-repeat rule) really
submitted a job (`qacct`-visible, external job ID confirmed), really ran
on a compute node, and reproduced the exact previously-verified counts
(203 ribosomal / 194 spliceosomal / 97 other) byte-for-byte -- not a
simulated success.

**Notebook rendering is a separate, local step, not part of this
Snakefile.** The notebooks live in the git repo on a laptop, not on
Hoffman2, and rendering them (pandas/matplotlib/seaborn via nbclient) is
cheap compared to the cluster-bound prediction stages -- there's little
to gain and real risk (a second conda env, git-on-a-cluster) in also
running that step on Hoffman2. After the `all` target finishes, run
`notebooks/sync_and_render.sh <ssh-host>` locally: it scp's the 7 tables
down and re-executes all 6 exploration notebooks in place.

**Out of scope for this workflow:** Deep-RBPPred's satmut array isn't
wrapped yet (it exists as a standalone `qsub` script,
`saturation_mutagenesis/submit_satmut_deeprbppred_array.sh`, but was
never merged into a combined output the way TANGO/protGPS/DeepLoc were,
so there's no stable final target to add as a rule). New rules can be
added to the Snakefile following the same pattern once that's resolved.

## Filtering down to RNA-binding proteins

Everything above runs genome-wide, on purpose (nothing about SpliceAI/
DeepLoc/protGPS is RBP-specific) -- but this project's actual scope is
RNA-binding proteins, so the six merged tables need a filtering step
before use, not to be read as-is. See `rbp_filter/README.md` for the RBP
gene list used, the row counts before and after, and the scripts.

## Saturation mutagenesis (RBPs only)

A separate, much larger scan -- every possible missense point mutant
of all 2,336 RBPs (~31.5M sequences), scored with DeepLoc + protGPS --
lives in `saturation_mutagenesis/README.md`. Uses a labmate's
pre-built mutant-sequence file as input rather than this repo's own
`mutant_sequence_construction/` scripts; see that README for why, and
for the GPU-scheduling issues (account-wide quota caps, a silent
task-failure mode, the queue-list fix) hit while running it at scale.

## A real failure worth documenting: don't trust a small-scale timing test

The first full-scale SpliceAI run (200-task CPU array job) was sized
using the per-variant rate observed in a 20-variant test
(~100-150ms/variant), giving an estimated ~47min/chunk and an
`h_rt=1:30:00` limit. Every task hit that wall-clock limit and was
killed (`execd enforced h_rt limit`, exit 137) without finishing — the
small test's rate didn't hold up on real (shared, contended) general-pool
nodes, and SpliceAI doesn't write its output VCF incrementally, so a
killed task loses all its computed work, not just the unfinished tail.
Confirmed via `qacct -j <job> -t <task>` showing `maxvmem` up to 17G
against a requested `h_data=4G`, and `failed: 44 : execd enforced h_rt
limit` on multiple sampled tasks. Fixed by rerunning with `h_rt=6:00:00`
and `h_data=24G`, and by giving each array task its own log file
(`-o .../log.$TASK_ID` — SGE does *not* auto-split a shared `-o` path
across array tasks, so all 200 tasks' stdout was interleaved into one
unreadable file the first time). Lesson: verify actual output row counts
after a "successful" (exit-0, all-tasks-finished) array job, not just
that every task exited without error — an `h_rt`-killed task can still
show as cleanly gone from `qstat` afterward.
