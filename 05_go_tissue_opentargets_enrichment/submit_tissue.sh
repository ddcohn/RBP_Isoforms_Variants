#!/bin/bash
#$ -N tissue_pull
#$ -cwd
#$ -o /u/home/d/ddcohn/_claude_qsub_tissue.log
#$ -j y
#$ -l h_data=4G,h_rt=6:00:00
#$ -pe shared 1

python3 /u/project/kappel/ddcohn/RBP/scripts/pull_tissue.py
