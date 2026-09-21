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
