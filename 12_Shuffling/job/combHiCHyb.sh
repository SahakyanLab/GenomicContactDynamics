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
################################################
####Job Details#################################
################################################

module load R/3.5.0-newgcc
module load gcc/4.9.2

R --vanilla < /t1-data/user/ltamon/DPhil/GenomicContactDynamics/8_ShuffleContactBins/B1.5_HicHybridisation_comb.R
