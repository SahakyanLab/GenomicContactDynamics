#!/bin/sh
##########################################################################
## A script t4mplate for submitting batch jobs.
## Please note that anything after the first two characters "#$" on a line
## will be treated as a SUN Grid Engine command.
##########################################################################
#$ -cwd
#$ -q batchq
#$ -l h_vmem=4G
#$ -pe dedicated 4
#$ -t 614-712
#$ -M ltamon
##$ -m eas
#########################################################################
## JOB DETAILS * JOB DETAILS * JOB DETAILS * JOB DETAILS * JOB DETAILS ##
#########################################################################
module load R/3.6.0-newgcc
module load gcc/4.9.2

Rscript --vanilla /t1-data/user/ltamon/DPhil/GenomicContactDynamics/20_ChromFeatAssoc/C1_associate/CptopCP3/script/assocCptopCP3${SGE_TASK_ID}.R
