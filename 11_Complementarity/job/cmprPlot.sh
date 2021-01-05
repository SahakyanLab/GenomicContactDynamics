#!/bin/sh
##########################################################################
## A script template for submitting batch jobs.
## Please note that anything after "#SBATCH" on a line will be treated as
## a SLURM command.
##########################################################################
#SBATCH -p batch
#SBATCH --mem-per-cpu=50G
#SBATCH -n 1
#SBATCH --cpus-per-task=1
#SBATCH --mail-user=ltamon
#SBATCH --mail-type=ALL
##########################################################################
module load R-base/4.0.1
module load R-cbrg/current

R --vanilla < /t1-data/user/ltamon/DPhil/GenomicContactDynamics/11_Constraints/D2_compare_plot.R
