#!/bin/sh
##########################################################################
## A script template for submitting batch jobs.
## Please note that anything after "#SBATCH" on a line will be treated as
## SLURM command.
##########################################################################
#SBATCH --mem-per-cpu=7G
#SBATCH -n 7
#SBATCH --cpus-per-task=7
#SBATCH --array=1-23
#########################################################################
## JOB DETAILS * JOB DETAILS * JOB DETAILS * JOB DETAILS * JOB DETAILS ##
#########################################################################
module load R-base/4.0.1
module load R-cbrg/current

R --vanilla < /t1-data/user/ltamon/DPhil/GenomicContactDynamics/21_Simulation/B0_contactMxGap/script/contactMxGap${SLURM_ARRAY_TASK_ID}.R

