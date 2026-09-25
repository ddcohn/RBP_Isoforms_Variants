#!/bin/bash -l
#$ -N merge_deeploc
#$ -cwd
#$ -o /u/home/d/ddcohn/_claude_qsub_merge_deeploc.log
#$ -j y
#$ -l h_data=48G,h_rt=3:00:00
#$ -t 1-2

module load miniforge/23.11.0
source activate spliceai
WHICH=(clinvar cmc)
python3 /u/home/d/ddcohn/_claude_merge_deeploc_deltas.py ${WHICH[$((SGE_TASK_ID - 1))]}
