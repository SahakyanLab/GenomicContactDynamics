#!/bin/sh
##########################################################################
## A script template for submitting batch jobs.
## Please note that anything after "#SBATCH" on a line will be treated as
## a SLURM command.
##########################################################################
#SBATCH -p batch
#SBATCH --mem-per-cpu=30G
#SBATCH -n 1
#SBATCH --mail-user=ltamon
#SBATCH --mail-type=ALL
##########################################################################
module load R-base/4.0.1
module load R-cbrg/current
module load gcc/9.3.0

R --vanilla < /t1-data/user/ltamon/DPhil/GenomicContactDynamics/19_Mutation_rates/A2_filter.R

