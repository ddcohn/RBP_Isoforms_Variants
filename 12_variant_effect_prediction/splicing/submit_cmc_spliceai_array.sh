#!/bin/bash -l
#$ -N spliceai_cmc
#$ -cwd
#$ -o /u/home/d/ddcohn/spliceai_task_logs/cmc.log.$TASK_ID
#$ -j y
#$ -l h_data=24G,h_rt=6:00:00
#$ -t 1-260
#$ -pe shared 1

module load miniforge/23.11.0
source activate spliceai

CHUNKDIR=/u/project/kappel/ddcohn/SpliceAI/cmc_run/chunks
OUTDIR=/u/project/kappel/ddcohn/SpliceAI/cmc_run/results
mkdir -p $OUTDIR

TASK=$(printf "%03d" $SGE_TASK_ID)
IN_VCF=$CHUNKDIR/chunk_${TASK}.vcf
OUT_VCF=$OUTDIR/chunk_${TASK}_out.vcf

if [ ! -f "$IN_VCF" ]; then
  echo "No input chunk for task $SGE_TASK_ID, skipping."
  exit 0
fi

START=$(date +%s)
spliceai -I $IN_VCF -O $OUT_VCF \
  -R /u/project/kappel/ddcohn/RNA-GPS/rnagps/reference/GRCh38.primary_assembly.genome.fa \
  -A grch38
END=$(date +%s)

echo "Task $SGE_TASK_ID done in $((END-START))s, $(wc -l < $OUT_VCF) output lines."
