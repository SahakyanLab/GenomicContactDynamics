#!/bin/sh
##########################################################################
## A script template for submitting batch jobs.
## Please note that anything after "#SBATCH" on a line will be treated as
## SLURM command.
##########################################################################
#SBATCH -p batch
#SBATCH --mem-per-cpu=4G
#SBATCH -n 5
#SBATCH --cpus-per-task=5
#SBATCH --array=721-732
#########################################################################
## JOB DETAILS * JOB DETAILS * JOB DETAILS * JOB DETAILS * JOB DETAILS ##
#########################################################################
module load R-base/4.0.1
module load R-cbrg/current
module load gcc/9.3.0

Rscript --vanilla /t1-data/user/ltamon/DPhil/GenomicContactDynamics/20_ChromFeatAssoc/C1_associate/CptopCP3/script/assocCptopCP3${SLURM_ARRAY_TASK_ID}.R
