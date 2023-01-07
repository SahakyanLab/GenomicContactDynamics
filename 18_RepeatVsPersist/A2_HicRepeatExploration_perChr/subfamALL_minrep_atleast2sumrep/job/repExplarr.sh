#!/bin/sh
##########################################################################
## A script template for submitting batch jobs.
## Please note that anything after the first two characters "#$" on a line
## will be treated as a SUN Grid Engine command.
##########################################################################
#SBATCH --mem=50G
#SBATCH -n 1
#SBATCH --array=21-22
#SBATCH --mail-user=ltamon
#SBATCH --mail-type=ALL
#########################################################################
## JOB DETAILS * JOB DETAILS * JOB DETAILS * JOB DETAILS * JOB DETAILS ##
#########################################################################
module load R-base/4.1.2
module load R-cbrg/current

R --vanilla < /project/sahakyanlab/ltamon/DPhil/GenomicContactDynamics/4_RepeatVsPersist/A2_HicRepeatExploration/subfamALL_minrep_atleast2sumrep/script/repExpl${SLURM_ARRAY_TASK_ID}.R

