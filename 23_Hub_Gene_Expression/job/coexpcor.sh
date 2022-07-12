#!/bin/sh
##########################################################################
## A script template for submitting batch jobs.
## Please note that anything after "#SBATCH" on a line will be treated as
## a SLURM command.
##########################################################################
#SBATCH --mem=60G
#SBATCH -n 5
#SBATCH --mail-user=ltamon
#SBATCH --mail-type=ALL
##########################################################################
module load R-base/4.1.0
module load R-cbrg/current

R --vanilla < /project/sahakyanlab/ltamon/SahakyanLab/CoreGenomeExplorer/F1_coexpression_pairCor.R
