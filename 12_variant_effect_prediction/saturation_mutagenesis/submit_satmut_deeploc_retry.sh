#!/bin/bash -l
#$ -N deeploc_satmut_retry
#$ -cwd
#$ -o /u/home/d/ddcohn/deeploc_task_logs/satmut_retry.log.$TASK_ID
#$ -j y
#$ -q gpu_RTX2080Ti.q,gpu_l40s.q@g18,gpu_l40s.q@g19,gpu_a100.q,gpu_h100.q,kappel_gpu.q
#$ -soft -q kappel_gpu.q,gpu_l40s.q,gpu_a100.q,gpu_h100.q
#$ -l h_data=24G,h_rt=8:00:00,cuda=1,require_gpu=1
#$ -t 1-5
#$ -pe shared 1

# Edit CHUNKS below to the task numbers that need redoing. Two distinct
# failure signatures have shown up in this pipeline, both diagnosed via
# qacct -j <job_id> -t <task_id>:
#   - Wall-time kill (contended/slow queue, chunk never finished in time):
#     exit_status=137, failed=44, "execd enforced h_rt limit".
#   - CUDA crash (illegal memory access inside DeepLoc2's model, not a
#     timeout): exit_status=1, failed=0, task finishes in well under a
#     minute, log shows "RuntimeError: CUDA error: an illegal memory
#     access was encountered". Seen across 10 different chunks in one
#     round; a plain resubmit (no special handling needed) succeeded for
#     all of them on the next attempt -- so far looks transient/hardware-
#     probabilistic, not deterministic on chunk content, but if it keeps
#     recurring on the exact same chunk across multiple resubmits that
#     assumption should be revisited.
# #$ -t above must match CHUNKS' length. This project's real backfill
# history (for reference, not meant to be re-run): round 1 caught chunks
# 0005/0006/0009/0010/0011 via wall-time kills; a later full audit found
# chunk 0006 had been silently dropped from an earlier resubmit (an
# actual scripting mistake, not a cluster issue) plus a fresh batch of 10
# CUDA-crash chunks and 1 more wall-time kill; all were cleared within
# two more rounds.
CHUNKS=(0005 0006 0009 0010 0011)
TASK=${CHUNKS[$((SGE_TASK_ID - 1))]}

module load miniforge/23.11.0
source activate deeploc2

CHUNKDIR=/u/project/kappel/ddcohn/protein_variant_effects/rbp_only/satmut_chunks
OUTDIR=/u/project/kappel/ddcohn/protein_variant_effects/rbp_only/satmut_deeploc_results

IN_FASTA=$CHUNKDIR/chunk_${TASK}.fasta
OUT_DIR=$OUTDIR/chunk_${TASK}

rm -rf $OUT_DIR
START=$(date +%s)
deeploc2 -f $IN_FASTA -o $OUT_DIR -d cuda
RC=$?
END=$(date +%s)

n_in=$(grep -c '^>' $IN_FASTA)
n_out=$(($(wc -l < $OUT_DIR/results_*.csv 2>/dev/null || echo 0) - 1))
echo "Retry task $SGE_TASK_ID (chunk $TASK): rc=$RC elapsed $((END-START))s input=$n_in output=$n_out"
if [ "$RC" -ne 0 ] || [ "$n_out" -lt "$n_in" ]; then
  exit 1
fi
