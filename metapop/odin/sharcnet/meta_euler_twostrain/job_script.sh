#!/bin/bash
#SBATCH --account=def-bolker
#SBATCH --job-name=meta_twostrain
#SBATCH --time=0-08:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/meta_twostrain_%j.out

module load r/4.5.0
task.run
