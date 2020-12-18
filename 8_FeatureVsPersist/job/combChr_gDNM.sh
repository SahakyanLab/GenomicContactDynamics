#!/bin/sh
##########################################################################
## A script template for submitting batch jobs.
## Please note that anything after "#SBATCH" on a line will be treated as
## SLURM command.
##########################################################################
#SBATCH -p batch
#SBATCH --mem-per-cpu=15G
#SBATCH -n 3
#SBATCH --cpus-per-task=3
#########################################################################
## JOB DETAILS * JOB DETAILS * JOB DETAILS * JOB DETAILS * JOB DETAILS ##
#########################################################################
module load R-base/4.0.1
module load R-cbrg/current
module load gcc/9.3.0

R --vanilla < /t1-data/user/ltamon/DPhil/GenomicContactDynamics/10_ChromatinFeatures/A3_combChr.R
