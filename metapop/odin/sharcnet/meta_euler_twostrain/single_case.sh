#!/bin/bash
# task.run calls: source ./single_case.sh $TABLE $i
# $TABLE ($1) = line content from table.dat = case ID (1..N_cases)
# $i ($2) = loop counter (same value here since table.dat is seq 1 14400)
# Must use a subshell for cd: task.run *sources* this script, so a bare
# cd would permanently change task.run's working directory and break
# subsequent cases.
module load r/4.5.0
export SLURM_ARRAY_TASK_ID=$1
(cd .. && Rscript euler_twostrain_run_array.R)
STATUS=$?
