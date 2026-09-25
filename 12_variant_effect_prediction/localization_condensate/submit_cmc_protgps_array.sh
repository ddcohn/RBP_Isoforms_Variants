#!/bin/bash -l
#$ -N protgps_cmc
#$ -cwd
#$ -o /u/home/d/ddcohn/protgps_task_logs/protgps_cmc.log.$TASK_ID
#$ -j y
#$ -q kappel_gpu.q
#$ -l h_data=24G,h_rt=8:00:00,cuda=1,require_gpu=1,highp
#$ -t 1-32
#$ -pe shared 1

module load miniforge/23.11.0
source activate protgps

TASK=$(printf "%03d" $SGE_TASK_ID)
python3 /u/home/d/ddcohn/_claude_protgps_predict_cmc_scale.py $TASK
