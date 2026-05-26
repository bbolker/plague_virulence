#!/bin/bash
#SBATCH --account=def-bolker
#SBATCH --job-name=euler_extinct
#SBATCH --array=1-280
#SBATCH --time=0-02:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=4
#SBATCH --output=logs/euler_%A_%a.out

# Submit from the metapop/odin directory so here::here() resolves correctly.
# Before first submission: mkdir -p logs outputs
# Check available R versions with: module spider r

module load r/4.5.0

Rscript euler_onepatch_onestrain_extinct_run_array.R
