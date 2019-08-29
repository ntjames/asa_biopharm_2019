#!/bin/bash
#SBATCH --job-name=cbmod_calib
#SBATCH --mail-user=nathan.t.james@vanderbilt.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=00:07:30
#SBATCH --mem=2G
#SBATCH --constraint=sandybridge|haswell
#SBATCH --array=1-5 
#SBATCH --partition=debug ## change to 'production' for actual run
#SBATCH --output=cbout/cbmod_sim_%a.out

module restore R_360
R --version

echo "SLURM_JOBID: " $SLURM_JOBID
echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
echo "SLURM_ARRAY_JOB_ID: " $SLURM_ARRAY_JOB_ID

## run this first to make simarray and compile Stan models
## Rscript --no-save 10_cb_calib_simarray.R
## Rscript --no-save 11_cbmod_calib_init.R

## run array job over sim ids
## there are 12150 total sim ids
Rscript --no-save 12_cb_calib_run.R $SLURM_ARRAY_TASK_ID
