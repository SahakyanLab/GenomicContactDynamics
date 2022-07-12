#!/bin/sh
##########################################################################
## A script template for submitting batch jobs.
## Please note that anything after the first two characters "#$" on a line
## will be treated as a SUN Grid Engine command.
##########################################################################
#SBATCH --mem=5G
#SBATCH -n 1
#SBATCH --array=1-23
#SBATCH --mail-user=ltamon
#SBATCH --mail-type=ALL
#########################################################################
## JOB DETAILS * JOB DETAILS * JOB DETAILS * JOB DETAILS * JOB DETAILS ##
#########################################################################
module load R-base/4.1.2
module load R-cbrg/current

R --vanilla < /project/sahakyanlab/ltamon/SahakyanLab/CoreGenomeExplorer/F1_coexpression_pairCor/script/pairCor${SLURM_ARRAY_TASK_ID}.R
