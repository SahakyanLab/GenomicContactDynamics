#!/bin/sh
##########################################################################
## A script template for submitting batch jobs.
## Please note that anything after the first two characters "#$" on a line
## will be treated as a SUN Grid Engine command.
##########################################################################
#$ -cwd
#$ -q batchq
#$ -l h_vmem=12G
#$ -pe dedicated 3
#$ -t 17

##$ -t 1-4 40G-3
##$ -t 23 30G-3
##$ -t 5-9 20G-3
##$ -t 10-17 12G-3
##$ -t 18-22 5G-3

#$ -M ltamon
#$ -m eas
#########################################################################
## JOB DETAILS * JOB DETAILS * JOB DETAILS * JOB DETAILS * JOB DETAILS ##
#########################################################################
module load R/3.6.0-newgcc
module load gcc/4.9.2

Rscript --vanilla /t1-data/user/ltamon/DPhil/GenomicContactDynamics/8_ShuffleContactBins/B1_getConstraints/script_align/getConstraints${SGE_TASK_ID}.R
