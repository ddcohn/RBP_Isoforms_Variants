#!/bin/bash
#$ -N combine_satmut
#$ -cwd
#$ -o /u/home/d/ddcohn/combine_satmut.log
#$ -j y
#$ -l h_data=32G,h_rt=6:00:00
#$ -pe shared 2

module load miniforge/23.11.0
python3 /u/home/d/ddcohn/_claude_combine_satmut_predictions.py
