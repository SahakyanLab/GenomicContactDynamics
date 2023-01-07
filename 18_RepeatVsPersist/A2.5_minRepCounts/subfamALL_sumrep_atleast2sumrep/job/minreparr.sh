#!/bin/sh
##########################################################################
## A script template for submitting batch jobs.
## Please note that anything after the first two characters "#$" on a line
## will be treated as a SUN Grid Engine command.
##########################################################################
#SBATCH --mem=100G
#SBATCH -n 1
#SBATCH --array=10-23 #1-2 200G #4-8 150G #9-12, 23 100G #13-22 50G
#SBATCH --mail-user=ltamon
#SBATCH --mail-type=ALL
#########################################################################
## JOB DETAILS * JOB DETAILS * JOB DETAILS * JOB DETAILS * JOB DETAILS ##
#########################################################################
module load R-base/4.2.0
module load R-cbrg/current

R --vanilla < /project/sahakyanlab/ltamon/SahakyanLab/GenomicContactDynamics/18_RepeatVsPersist/A2.5_minRepCounts/subfamALL_sumrep_atleast2sumrep/script/minrep${SLURM_ARRAY_TASK_ID}.R
