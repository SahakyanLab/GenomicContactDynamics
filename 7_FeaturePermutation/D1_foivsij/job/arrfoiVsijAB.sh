#!/bin/sh
##########################################################################
## A script t4mplate for submitting batch jobs.
## Please note that anything after the first two characters "#$" on a line
## will be treated as a SLURM command.
##########################################################################
#SBATCH --mem=8G
#SBATCH -n 1
#SBATCH --array=125-189
##SBATCH --mail-user=ltamon
##SBATCH --mail-type=ALL
#########################################################################
## JOB DETAILS * JOB DETAILS * JOB DETAILS * JOB DETAILS * JOB DETAILS ##
#########################################################################
module load R-base/4.1.0
module load R-cbrg/current

R --vanilla < /stopgap/sahakyanlab/ltamon/DPhil/GenomicContactDynamics/20_ChromFeatAssoc/D1_foivsij/script/any/foivsij${SLURM_ARRAY_TASK_ID}.R
