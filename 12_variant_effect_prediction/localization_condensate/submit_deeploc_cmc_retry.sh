#!/bin/bash -l
#$ -N deeploc_cmc_retry
#$ -cwd
#$ -o /u/home/d/ddcohn/deeploc_task_logs/deeploc_cmc_retry.log.$TASK_ID
#$ -j y
#$ -q gpu_l40s.q
#$ -l h_data=24G,h_rt=8:00:00,cuda=1,require_gpu=1
#$ -t 1-2
#$ -pe shared 1

module load miniforge/23.11.0
source activate deeploc2

CHUNKS=(009 017)
TASK=${CHUNKS[$((SGE_TASK_ID - 1))]}

IN_FASTA=/u/project/kappel/ddcohn/protein_variant_effects/cmc_deeploc_chunks/chunk_${TASK}.fasta
OUT_DIR=/u/project/kappel/ddcohn/protein_variant_effects/cmc_deeploc_results_backup/chunk_${TASK}

nvidia-smi --query-gpu=index,name,memory.total --format=csv
START=$(date +%s)
deeploc2 -f $IN_FASTA -o $OUT_DIR -d cuda
RC=$?
END=$(date +%s)

echo "Retry task $SGE_TASK_ID (chunk $TASK) exit code $RC, elapsed $((END-START))s."
exit $RC
