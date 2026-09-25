#!/bin/bash -l
#$ -N deeploc_cmc
#$ -cwd
#$ -o /u/home/d/ddcohn/deeploc_task_logs/deeploc_cmc.log.$TASK_ID
#$ -j y
#$ -l h_data=24G,h_rt=20:00:00,cuda=1,require_gpu=1
#$ -t 1-32
#$ -pe shared 1

module load miniforge/23.11.0
source activate deeploc2

CHUNKDIR=/u/project/kappel/ddcohn/protein_variant_effects/cmc_deeploc_chunks
OUTDIR=/u/project/kappel/ddcohn/protein_variant_effects/cmc_deeploc_results
mkdir -p $OUTDIR

TASK=$(printf "%03d" $SGE_TASK_ID)
IN_FASTA=$CHUNKDIR/chunk_${TASK}.fasta
OUT_DIR=$OUTDIR/chunk_${TASK}

if [ ! -f "$IN_FASTA" ]; then
  echo "No input chunk for task $SGE_TASK_ID, skipping."
  exit 0
fi

nvidia-smi --query-gpu=index,memory.used,memory.total --format=csv
START=$(date +%s)
deeploc2 -f $IN_FASTA -o $OUT_DIR -d cuda
RC=$?
END=$(date +%s)

echo "Task $SGE_TASK_ID exit code $RC, elapsed $((END-START))s."
exit $RC
