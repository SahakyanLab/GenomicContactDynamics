#!/bin/sh
##########################################################################
## A script template for submitting batch jobs.
## Please note that anything after the first two characters "#$" on a line
## will be treated as a SUN Grid Engine command.
##########################################################################
#$ -cwd
#$ -q batchq
#$ -l h_vmem=15G
#$ -pe dedicated 3

##$ -t 18-22 5G-3
##$ -t 2 50G-3
##$ -t 10-17 15G-3
##$ -t 4-9 20G-3
##$ -t 23 20G-3
##$ -t 3-4 30G-3
##$ -t 1 40G-3
#$ -t 17

#$ -M ltamon
#$ -m eas
#########################################################################
## JOB DETAILS * JOB DETAILS * JOB DETAILS * JOB DETAILS * JOB DETAILS ##
#########################################################################
module load R/3.6.0-newgcc
module load gcc/4.9.2

Rscript --vanilla /t1-data/user/ltamon/DPhil/GenomicContactDynamics/11_Constraints/B2_getConstraints/script/B2_getConstraints${SGE_TASK_ID}.R
