#!/bin/sh
##########################################################################
## A script template for submitting batch jobs.
## Please note that anything after the first two characters "#$" on a line
## will be treated as a SUN Grid Engine command.
##########################################################################
#$ -cwd
#$ -q batchq
#$ -M ltamon
#$ -m eas
#########################################################################
## JOB DETAILS * JOB DETAILS * JOB DETAILS * JOB DETAILS * JOB DETAILS ##
#########################################################################
module load R/3.6.0-newgcc
module load gcc/4.9.2

# --vanilla combines --no-save, --no-restore, --no-site-file, --no-init-file and --no-environ
Rscript --vanilla /t1-data/user/ltamon/DPhil/GenomicContactDynamics/11_Constraints/B2.5_fixConstraints.R
