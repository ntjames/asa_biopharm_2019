#!/bin/bash
#SBATCH --job-name=sbcmodsim
#SBATCH --mail-user=nathan.t.james@vanderbilt.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:05:00
#SBATCH --mem=1G
#SBATCH --constraint=sandybridge|haswell
#SBATCH --array=1-1000
#SBATCH --output=sbc_out/sbcmod_sim_%a.out

module restore R_basic
R --version 

echo "SLURM_JOBID: " $SLURM_JOBID
echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
echo "SLURM_ARRAY_JOB_ID: " $SLURM_ARRAY_JOB_ID

## run this first to compile STAN models
## Rscript --no-save sbc_mod_init.R

## run array job over sim ids
##Rscript --no-save sbc_mod_accre.R $SLURM_ARRAY_TASK_ID
Rscript sbc_mod_accre.R $SLURM_ARRAY_TASK_ID
