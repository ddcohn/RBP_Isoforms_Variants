#!/bin/bash
#$ -N cosmic_table
#$ -cwd
#$ -o /u/home/d/ddcohn/_claude_qsub_cosmic_table.log
#$ -j y
#$ -l h_data=8G,h_rt=2:00:00
#$ -pe shared 1

python3 /u/home/d/ddcohn/_claude_build_cosmic_table.py
