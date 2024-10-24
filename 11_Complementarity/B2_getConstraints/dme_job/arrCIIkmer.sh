#!/bin/sh
##########################################################################
## A script template for submitting batch jobs.
## Please note that anything after "#SBATCH" on a line will be treated as
## a SLURM command.
##########################################################################
#SBATCH --mem=10G
#SBATCH -n 3
#SBATCH --array=1-8
##SBATCH --mail-user=ltamon
##SBATCH --mail-type=ALL
##########################################################################
module load R-base/4.1.0
module load R-cbrg/current

R --vanilla < /project/sahakyanlab/ltamon/DPhil/GenomicContactDynamics/11_Constraints/B2_getConstraints/dme_script/chr${SLURM_ARRAY_TASK_ID}typekmer.R
