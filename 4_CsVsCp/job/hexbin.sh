#!/bin/sh
##########################################################################
## A script template for submitting batch jobs.
## Please note that anything after "#SBATCH" on a line will be treated as
## a SLURM command.
##########################################################################
#SBATCH -n 1
#SBATCH --mail-user=ltamon
#SBATCH --mail-type=ALL
##########################################################################
module load R-base/4.0.1
module load R-cbrg/current

R --vanilla < /project/sahakyanlab/ltamon/DPhil/GenomicContactDynamics/4_CsVsCp/A2_hexbinplot.R
