#!/bin/bash
#SBATCH --job-name=bbmod_calib
#SBATCH --mail-user=nathan.t.james@vanderbilt.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=00:05:30
#SBATCH --mem=2G
#SBATCH --constraint=sandybridge|haswell
#SBATCH --array=31-500
#SBATCH --output=bbout/bbmod_sim_%a.out

module restore R_360
R --version

echo "SLURM_JOBID: " $SLURM_JOBID
echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
echo "SLURM_ARRAY_JOB_ID: " $SLURM_ARRAY_JOB_ID

## run this first to compile Stan models
## Rscript --no-save 0_bbmod_calib_init.R

## run array job over sim ids
## there are 12150 total sim ids
Rscript --no-save 1_bb_calib_run.R $SLURM_ARRAY_TASK_ID
