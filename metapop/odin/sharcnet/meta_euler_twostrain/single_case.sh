#!/bin/bash
# task.run calls: source ./single_case.sh $TABLE $i
# $1 ($TABLE) = path to table.dat
# $2 ($i)     = case index (line number in table.dat)
#
# Must use a subshell for cd: task.run *sources* this script, so a bare
# cd would permanently change task.run's working directory and break
# subsequent cases.
TABLE1=$1
i1=$2

LINE=$(sed -n ${i1}p "$TABLE1")
ID=$(echo "$LINE" | cut -d" " -f1)

# Capture metajob ID before we overwrite SLURM_ARRAY_TASK_ID below:
METAJOB_ID=${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}

module load r/4.5.0
(cd .. && export SLURM_ARRAY_TASK_ID=$ID && Rscript euler_twostrain_run_array.R)
STATUS=$?

echo "$ID $STATUS" >> STATUSES/status.$METAJOB_ID
