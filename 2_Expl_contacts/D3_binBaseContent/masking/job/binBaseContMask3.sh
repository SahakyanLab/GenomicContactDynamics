#!/bin/sh
##########################################################################
## A script template for submitting batch jobs.
## Please note that anything after "#SBATCH" on a line will be treated as
## a SLURM command.
##########################################################################
#SBATCH --mem=5G
#SBATCH -n 1
#SBATCH --array=1-23
##SBATCH --mail-user=ltamon
##SBATCH --mail-type=ALL
##########################################################################
module load R-base/4.0.1
module load R-cbrg/current

Rscript --vanilla /t1-data/user/ltamon/DPhil/GenomicContactDynamics/2_Expl_contacts/D3_binBaseContent/masking/script/chr${SLURM_ARRAY_TASK_ID}mask3.R
