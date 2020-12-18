#!/bin/sh
##########################################################################
## A script template for submitting batch jobs.
## Please note that anything after "#SBATCH" on a line will be treated as
## a SLURM command.
##########################################################################
#SBATCH -p batch
#SBATCH --mem-per-cpu=10G
#SBATCH -n 1
#SBATCH --cpus-per-task=1
#SBATCH --mail-user=ltamon
#SBATCH --mail-type=ALL
##########################################################################
module load R-base/4.0.1
module load R-cbrg/current

Rscript --vanilla /t1-data/user/ltamon/DPhil/GenomicContactDynamics/6_Chrom3D/C4_ContactRadDist_plot.R
