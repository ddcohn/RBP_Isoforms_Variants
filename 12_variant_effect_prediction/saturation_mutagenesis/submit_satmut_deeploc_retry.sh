#!/bin/bash -l
#$ -N deeploc_satmut_retry
#$ -cwd
#$ -o /u/home/d/ddcohn/deeploc_task_logs/satmut_retry.log.$TASK_ID
#$ -j y
#$ -q gpu_RTX2080Ti.q,gpu_l40s.q,gpu_a100.q,gpu_h100.q
#$ -soft -q gpu_l40s.q,gpu_a100.q,gpu_h100.q
#$ -l h_data=24G,h_rt=8:00:00,cuda=1,require_gpu=1
#$ -t 1-5
#$ -pe shared 1

# Edit CHUNKS below to the task numbers that need redoing (see README:
# "GPU scheduling failures" for how to find them -- qacct exit_status=137
# with failed=44, "execd enforced h_rt limit", is the signature of a
# chunk landing on a contended/slow queue and getting killed before
# finishing). #$ -t above must match CHUNKS' length.
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
