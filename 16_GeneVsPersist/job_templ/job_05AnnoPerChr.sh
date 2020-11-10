#!/bin/sh
##########################################################################
## A script template for submitting batch jobs. To submit a batch job,
## please type
##
##    qsub myprog.sh
##
## Please note that anything after the first two characters "#$" on a line
## will be treated as a SUN Grid Engine command.
##########################################################################

# The following to run programs in the current working directory
#$ -cwd

# Specify a queue
#$ -q batchq


# The following two lines will send an email notification when
# job is Ended/Aborted/Suspended - Please replace "UserName" with your username.
#
#$ -M ltamon
#$ -m eas

Rscript /t1-data/user/ltamon/DPhil/GenomicContactDynamics/3_AnnotationVsPersist/A3_anno_perChr_05.R
