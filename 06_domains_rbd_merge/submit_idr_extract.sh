#!/bin/bash
#$ -N idr_extract
#$ -cwd
#$ -o /u/home/d/ddcohn/_claude_qsub_idr_extract.log
#$ -j y
#$ -l h_data=16G,h_rt=1:00:00
#$ -pe shared 1

python3 /u/home/d/ddcohn/_claude_extract_idr.py
