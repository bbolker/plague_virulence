#!/bin/bash
#SBATCH --account=def-bolker
#SBATCH --job-name=twostrain_extinct
#SBATCH --array=1-1000
#SBATCH --time=0-01:00:00
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/twostrain_%A_%a.out

# Submit from the metapop/odin directory so here::here() resolves correctly.
# Before first submission: mkdir -p logs outputs
# Check available R versions with: module spider r

module load r/4.5.0

Rscript discrete_onepatch_twostrain_extinct_run_array.R
