#!/bin/sh
##########################################################################
## A script template for submitting batch jobs.
## Please note that anything after "#SBATCH" on a line will be treated as
## SLURM command.
##########################################################################
#SBATCH -p batch
#SBATCH --mem-per-cpu=10G
#SBATCH -n 4
#SBATCH --cpus-per-task=4
#SBATCH --array=1-23
#########################################################################
## JOB DETAILS * JOB DETAILS * JOB DETAILS * JOB DETAILS * JOB DETAILS ##
#########################################################################
module load R-base/4.0.1
module load R-cbrg/current
module load gcc/9.3.0

Rscript --vanilla /t1-data/user/ltamon/DPhil/GenomicContactDynamics/11_Complementarity/B2_getConstraints/script_kmer3/B2_getConstraints${SLURM_ARRAY_TASK_ID}.R
