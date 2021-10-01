#!/bin/sh
##########################################################################
## A script template for submitting batch jobs.
## Please note that anything after "#SBATCH" on a line will be treated as
## a SLURM command.
##########################################################################
#SBATCH --mem=30G
#SBATCH -n 1
#SBATCH --mail-user=ltamon
#SBATCH --mail-type=ALL
##########################################################################
module load R-base/4.1.0
module load R-cbrg/current

R --vanilla < /stopgap/sahakyanlab/ltamon/DPhil/GenomicContactDynamics/pending/11_Constraints/D2_compare_plot_scatter.R
