#!/bin/sh
##########################################################################
## A script template for submitting batch jobs.
## Please note that anything after "#SBATCH" on a line will be treated as
## SLURM command.
##########################################################################
#SBATCH --mem=15G
#SBATCH -n 3
#SBATCH --array=733
#SBATCH --mail-user=ltamon
#SBATCH --mail-type=ALL
#########################################################################
## JOB DETAILS * JOB DETAILS * JOB DETAILS * JOB DETAILS * JOB DETAILS ##
#########################################################################
module load R-base/4.1.0
module load R-cbrg/current

Rscript --vanilla /project/sahakyanlab/ltamon/DPhil/GenomicContactDynamics/20_ChromFeatAssoc/C1_associate/CpAllCs1perc/script/assocCpAllCs1perc${SLURM_ARRAY_TASK_ID}.R
