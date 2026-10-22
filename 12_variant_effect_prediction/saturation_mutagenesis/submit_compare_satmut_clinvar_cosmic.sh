#!/bin/bash
#$ -N compare_satmut_clinvar_cosmic
#$ -cwd
#$ -o /u/home/d/ddcohn/compare_satmut_clinvar_cosmic.log
#$ -j y
#$ -l h_data=16G,h_rt=3:00:00
#$ -pe shared 1

module load miniforge/23.11.0
python3 /u/project/kappel/ddcohn/RBP_pipeline_scripts/saturation_mutagenesis/compare_satmut_clinvar_cosmic.py
