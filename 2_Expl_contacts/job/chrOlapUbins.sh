#!/bin/sh
##########################################################################
## A script template for submitting batch jobs.
## Please note that anything after the first two characters "#$" on a line
## will be treated as a SUN Grid Engine command.
##########################################################################
#$ -cwd
#$ -l h_vmem=4G
#$ -pe dedicated 3
#$ -M ltamon
#$ -m eas
#########################################################################
## JOB DETAILS * JOB DETAILS * JOB DETAILS * JOB DETAILS * JOB DETAILS ##
#########################################################################
module load R/3.6.0-newgcc
module load gcc/4.9.2

Rscript /t1-data/user/ltamon/DPhil/GenomicContactDynamics/2_HiC_Human21_Expl_ext/B4_ubinsOlapChr.R
