#!/bin/sh
##########################################################################
## A script template for submitting batch jobs.
## Please note that anything after the first two characters "#$" on a line
## will be treated as a SUN Grid Engine command.
##########################################################################

#$ -cwd
#$ -q batchq

#$ -l h_vmem=25G
#$ -pe dedicated 5

#$ -t 1-6

##$ -M ltamon
##$ -m eas

#########################################################################
## JOB DETAILS * JOB DETAILS * JOB DETAILS * JOB DETAILS * JOB DETAILS ##
#########################################################################
module load homer/liezel
module load R/3.5.0-newgcc
module load gcc/4.9.2

Rscript /t1-data/user/ltamon/DPhil/GenomicContactDynamics/15_Motif/A3_homer/A3_homer_bgrCp1_${SGE_TASK_ID}.R
