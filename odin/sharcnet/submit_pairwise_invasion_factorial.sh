#!/bin/bash
#SBATCH --account=def-bolker
#SBATCH --job-name=pairwise_factorial
#SBATCH --array=1-1323
#SBATCH --time=0-01:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/pairwise_factorial_%A_%a.out

# Submit from the odin/sharcnet directory, through submit_chunked.sh so that
# at most 50 tasks are queued or running at once (its --array overrides the
# one above):
#   bash submit_chunked.sh submit_pairwise_invasion_factorial.sh 1 1323 10 50 120
# Before first submission: mkdir -p logs outputs

module load r/4.5.0

export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

Rscript pairwise_invasion_factorial_run_array.R
