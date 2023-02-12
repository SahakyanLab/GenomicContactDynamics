#!/bin/sh
##########################################################################
## A script template for submitting batch jobs.
## Please note that anything after "#SBATCH" on a line will be treated as
## a SLURM command.
##########################################################################
#SBATCH --mem=20G
#SBATCH -n 1
#SBATCH --array=1-1400
##SBATCH --mail-user=ltamon
##SBATCH --mail-type=ALL
##########################################################################
module load R-base/4.2.0
module load R-cbrg/current

R --vanilla < /project/sahakyanlab/ltamon/SahakyanLab/GenomicContactDynamics/19_MutationRatesVsPersist/A2_plotdataVsCp_median/script/calcCALCREPLACEcomb${SLURM_ARRAY_TASK_ID}.R
