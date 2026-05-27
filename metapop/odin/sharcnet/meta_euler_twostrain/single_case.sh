#!/bin/bash
# $ID is the META-Farm case number (1-indexed line number from table.dat),
# which maps directly to SLURM_ARRAY_TASK_ID expected by the R script.
module load r/4.5.0
export SLURM_ARRAY_TASK_ID=$ID
cd ..
Rscript euler_twostrain_run_array.R
STATUS=$?
