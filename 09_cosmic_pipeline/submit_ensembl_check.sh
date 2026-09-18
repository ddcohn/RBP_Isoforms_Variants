#!/bin/bash
#$ -N ensembl_id_check
#$ -cwd
#$ -o /u/home/d/ddcohn/_claude_qsub_ensembl_check.log
#$ -j y
#$ -l h_data=4G,h_rt=3:00:00
#$ -pe shared 1

python3 /u/home/d/ddcohn/_claude_check_ensembl_id_full.py
