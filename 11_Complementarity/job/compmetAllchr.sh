#!/bin/sh
##########################################################################
## A script template for submitting batch jobs.
## Please note that anything after "#SBATCH" on a line will be treated as
## a SLURM command.
##########################################################################
#SBATCH --mem=40G
#SBATCH -n 3
#SBATCH --mail-user=ltamon
#SBATCH --mail-type=ALL
##########################################################################
module load R-base/4.1.0
module load R-cbrg/current

R --vanilla < /stopgap/sahakyanlab/ltamon/DPhil/GenomicContactDynamics/pending/11_Constraints/D4_compare_metric_allChr.R
