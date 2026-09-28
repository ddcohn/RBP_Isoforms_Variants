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

**One-time setup on Hoffman2:**
1. The `spliceai`, `deeploc2`, and `protgps` conda envs already need to
   exist (unchanged from before -- this workflow doesn't touch them, on
   purpose: they're GPU/version-pinned and already proven).
2. Create a `snakemake` env: `conda create -n snakemake -c bioconda -c
   conda-forge snakemake=7.32 -y`.
3. Place the raw input files at the paths in `workflow/config.yaml`
   (`ClinVar_variant_summary_complete.csv`,
   `CancerMutationCensus_AllData_v104_GRCh37.tsv.gz`) and the reference
   genome FASTA.

**Run it:**
```
module load miniforge/23.11.0 && source activate snakemake
cd 12_variant_effect_prediction/workflow
snakemake --profile profiles/hoffman2 all
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
`snakemake -n --profile profiles/hoffman2 all` dry-runs the whole DAG
without submitting anything (config.yaml is loaded automatically by the
Snakefile itself -- don't also pass `--configfile` on the command line,
Snakemake's `--configfile` flag greedily consumes the next argument too
and will misparse the target).
without submitting anything.

Target `all` produces the 6 merged delta/score tables (SpliceAI/DeepLoc/
protGPS x ClinVar/CMC) plus `clinvar_metadata_lookup.tsv`, all in
`/u/project/kappel/ddcohn/protein_variant_effects/`.

**Notebook rendering is a separate, local step, not part of this
Snakefile.** The notebooks live in the git repo on a laptop, not on
Hoffman2, and rendering them (pandas/matplotlib/seaborn via nbclient) is
cheap compared to the cluster-bound prediction stages -- there's little
to gain and real risk (a second conda env, git-on-a-cluster) in also
running that step on Hoffman2. After the `all` target finishes, run
`notebooks/sync_and_render.sh <ssh-host>` locally: it scp's the 7 tables
down and re-executes all 6 exploration notebooks in place.

**Out of scope for this workflow:** PTM (MusiteDeep) and protein
stability (ThermoMPNN) aren't started yet (see above), so there's
nothing to wrap for them. New rules can be added to the Snakefile
following the same pattern once those tools exist.

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
