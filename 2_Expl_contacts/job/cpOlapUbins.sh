#!/bin/sh
##########################################################################
## A script template for submitting batch jobs.
## Please note that anything after the first two characters "#$" on a line
## will be treated as a SUN Grid Engine command.
##########################################################################

#$ -cwd
#$ -M ltamon
#$ -m eas

#########################################################################
################################################
####Job Details#################################
################################################

module load R/3.5.0-newgcc
module load gcc/4.9.2

Rscript /t1-data/user/ltamon/DPhil/GenomicContactDynamics/2_HiC_Human21_Expl_ext/B3_ubinsOlapAcrossCps.R
