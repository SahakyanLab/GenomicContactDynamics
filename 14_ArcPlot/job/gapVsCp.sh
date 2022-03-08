#!/bin/sh
##########################################################################
## A script template for submitting batch jobs.
## Please note that anything after "#SBATCH" on a line will be treated as
## a SLURM command.
##########################################################################
#SBATCH --mem=50G
#SBATCH -n 1
#SBATCH --mail-user=ltamon
#SBATCH --mail-type=ALL
##########################################################################
module load R-base/4.1.2
module load R-cbrg/current

Rscript --vanilla /stopgap/sahakyanlab/ltamon/DPhil/GenomicContactDynamics/14_ArcPlot/A1_gapVsCp.R

