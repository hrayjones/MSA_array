#!/usr/bin/env bash
#SBATCH --mem=12GB

cd ~/MSA_array/scripts
Rscript 01_get_CG_proximal_SNPs.R

