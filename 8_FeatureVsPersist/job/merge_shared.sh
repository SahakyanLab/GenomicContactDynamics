#!/bin/sh
##########################################################################
## A script template for submitting batch jobs.
## Please note that anything after "#SBATCH" on a line will be treated as
## SLURM command.
##########################################################################
#SBATCH --mail-user=ltamon
#SBATCH --mail-type=ALL
#########################################################################
## JOB DETAILS * JOB DETAILS * JOB DETAILS * JOB DETAILS * JOB DETAILS ##
#########################################################################
module load R-base/4.1.0
module load R-cbrg/current

R --vanilla < /project/sahakyanlab/ltamon/DPhil/GenomicContactDynamics/10_ChromatinFeatures/A2_mergeFETA.R
