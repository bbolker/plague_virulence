#!/bin/bash
#SBATCH --account=def-bolker
#SBATCH --job-name=pairwise_factorial_mini
#SBATCH --array=1-32
#SBATCH --time=0-00:30:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/pairwise_factorial_mini_%A_%a.out

# Submit from the odin/sharcnet directory.
# Before first submission: mkdir -p logs outputs

module load r/4.5.0

export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

Rscript pairwise_invasion_factorial_run_array.R --mini
