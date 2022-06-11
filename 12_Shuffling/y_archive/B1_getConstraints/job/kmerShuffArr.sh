#!/bin/sh
##########################################################################
## A script template for submitting batch jobs.
## Please note that anything after the first two characters "#$" on a line
## will be treated as a SUN Grid Engine command.
##########################################################################
#$ -cwd
#$ -q batchq

#$ -l h_vmem=5G
#$ -pe dedicated 4

#$ -t 3-23 
##$ -t 1-2 10G-3

#$ -M ltamon
#$ -m eas
#########################################################################
## JOB DETAILS * JOB DETAILS * JOB DETAILS * JOB DETAILS * JOB DETAILS ##
#########################################################################
module load R/3.6.0-newgcc
module load gcc/4.9.2

Rscript --vanilla /t1-data/user/ltamon/DPhil/GenomicContactDynamics/8_ShuffleContactBins/B1_getConstraints/script_kmer/getConstraints${SGE_TASK_ID}.R
