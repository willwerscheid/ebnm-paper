#!/bin/bash
#SBATCH --job-name=ebnm-paper
#SBATCH --output=code_output_appendix.txt
#SBATCH --error=code_output_appendix.err
#SBATCH --time=36:00:00
#SBATCH --partition=broadwl
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G

module load R/3.6.1
Rscript code.R appendix-only