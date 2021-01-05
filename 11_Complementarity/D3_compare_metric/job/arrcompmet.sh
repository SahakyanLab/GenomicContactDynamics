#!/bin/sh
##########################################################################
## A script t4mplate for submitting batch jobs.
## Please note that anything after the first two characters "#$" on a line
## will be treated as a SLURM command.
##########################################################################
#SBATCH -p batch
#SBATCH --mem-per-cpu=12G
#SBATCH -n 1
#SBATCH --cpus-per-task=1
#SBATCH --array=1-6
#########################################################################
## JOB DETAILS * JOB DETAILS * JOB DETAILS * JOB DETAILS * JOB DETAILS ##
#########################################################################
module load R-base/4.0.1
module load R-cbrg/current

Rscript --vanilla /t1-data/user/ltamon/DPhil/GenomicContactDynamics/11_Constraints/D3_compare_metric/script/compmet${SLURM_ARRAY_TASK_ID}.R
