#!/bin/bash
#SBATCH -t 10
#SBATCH --array 1-2

echo "I am job $SLURM_JOBID, task $SLURM_ARRAY_TASK_ID"
