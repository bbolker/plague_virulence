#!/bin/bash
#SBATCH --account=def-bolker
#SBATCH --job-name=euler_twostrain_ex
#SBATCH --array=1-6
#SBATCH --time=0-02:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/euler_twostrain_examples_%A_%a.out

# Submit from the odin/sharcnet directory.
# Before first submission: mkdir -p logs outputs

module load r/4.5.0

Rscript euler_twostrain_examples_run_array.R
