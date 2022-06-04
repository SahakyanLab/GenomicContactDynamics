#!/bin/sh
##########################################################################
## A script template for submitting batch jobs.
## Please note that anything after "#SBATCH" on a line will be treated as
## a SLURM command.
##########################################################################
#SBATCH --mem=60G
#SBATCH -n 3
#SBATCH --array=101
#SBATCH --mail-user=ltamon
#SBATCH --mail-type=ALL
##########################################################################
module load R-base/4.1.2
module load R-cbrg/current

Rscript --vanilla /project/sahakyanlab/ltamon/DPhil/GenomicContactDynamics/12_MaskingFeatures/A1_masking/script/A1_maskedKmerCounts{SLURM_ARRAY_TASK_ID}.R
