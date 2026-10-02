#!/bin/bash -l
#$ -N deeploc_satmut
#$ -cwd
#$ -o /u/home/d/ddcohn/deeploc_task_logs/satmut.log.$TASK_ID
#$ -j y
#$ -q gpu_RTX2080Ti.q,gpu_l40s.q,gpu_a100.q,gpu_h100.q
#$ -soft -q gpu_l40s.q,gpu_a100.q,gpu_h100.q
#$ -l h_data=24G,h_rt=8:00:00,cuda=1,require_gpu=1
#$ -t 1-226
#$ -tc 15
#$ -pe shared 1

module load miniforge/23.11.0
source activate deeploc2

CHUNKDIR=/u/project/kappel/ddcohn/protein_variant_effects/rbp_only/satmut_chunks
OUTDIR=/u/project/kappel/ddcohn/protein_variant_effects/rbp_only/satmut_deeploc_results
mkdir -p $OUTDIR

TASK=$(printf "%04d" $SGE_TASK_ID)
IN_FASTA=$CHUNKDIR/chunk_${TASK}.fasta
OUT_DIR=$OUTDIR/chunk_${TASK}

if [ ! -f "$IN_FASTA" ]; then
  echo "No input chunk for task $SGE_TASK_ID, skipping."
  exit 0
fi

rm -rf $OUT_DIR
START=$(date +%s)
deeploc2 -f $IN_FASTA -o $OUT_DIR -d cuda
RC=$?
END=$(date +%s)

n_in=$(grep -c '^>' $IN_FASTA)
n_out=$(($(wc -l < $OUT_DIR/results_*.csv 2>/dev/null || echo 0) - 1))
echo "Task $SGE_TASK_ID (chunk $TASK): rc=$RC elapsed $((END-START))s input=$n_in output=$n_out"
if [ "$RC" -ne 0 ] || [ "$n_out" -lt "$n_in" ]; then
  exit 1
fi
