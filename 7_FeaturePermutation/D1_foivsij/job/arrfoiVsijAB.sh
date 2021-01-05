#!/bin/sh
##########################################################################
## A script t4mplate for submitting batch jobs.
## Please note that anything after the first two characters "#$" on a line
## will be treated as a SLURM command.
##########################################################################
#SBATCH -p batch
#SBATCH --mem-per-cpu=8G
#SBATCH -n 1
#SBATCH --cpus-per-task=1
#SBATCH --array=1-34
##SBATCH --mail-user=ltamon
##SBATCH --mail-type=ALL
#########################################################################
## JOB DETAILS * JOB DETAILS * JOB DETAILS * JOB DETAILS * JOB DETAILS ##
#########################################################################
module load R-base/4.0.1
module load R-cbrg/current
module load gcc/9.3.0

Rscript --vanilla /t1-data/user/ltamon/DPhil/GenomicContactDynamics/20_ChromFeatAssoc/D1_foivsij/script/any/foivsij${SLURM_ARRAY_TASK_ID}.R
