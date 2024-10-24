#!/bin/sh
##########################################################################
## A script template for submitting batch jobs.
## Please note that anything after "#SBATCH" on a line will be treated as
## a SLURM command.
##########################################################################
#SBATCH --mem=30G
#SBATCH -n 4
#SBATCH --mail-user=ltamon
#SBATCH --mail-type=ALL
##########################################################################
module load R-base/4.2.0
module load R-cbrg/current

R --vanilla < /project/sahakyanlab/ltamon/SahakyanLab/GenomicContactDynamics/19_MutationRatesVsPersist/A3.5_mut_contact_Cp_ACTUALplot/script/A3.5_mut_contact_Cp_ACTUALplot_icgc.R
