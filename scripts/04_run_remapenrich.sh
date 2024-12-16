#!/usr/bin/env bash
#SBATCH --mem=12GB

# If using /bin/bash above, do the following:
source ~/.bashrc
source ~/miniconda3/etc/profile.d/conda.sh # because the pipeline otherwise cannot activate the env

conda activate remapenrich

cd ~/MSA_array/scripts
Rscript 04_remapenrich.R
