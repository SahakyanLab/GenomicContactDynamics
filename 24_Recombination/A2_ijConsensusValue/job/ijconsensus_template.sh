#!/bin/sh
##########################################################################
## A script template for submitting batch jobs.
## Please note that anything after "#SBATCH" on a line will be treated as
## a SLURM command.
##########################################################################
#SBATCH --mem=15G
#SBATCH -n 3
#SBATCH --array=1-23
#SBATCH --mail-user=ltamon
#SBATCH --mail-type=ALL
##########################################################################
module load R-base/4.1.2
module load R-cbrg/current

R --vanilla < /project/sahakyanlab/ltamon/SahakyanLab/GenomicContactDynamics/24_Recombination/A2_ijConsensusValue/script/chr${SLURM_ARRAY_TASK_ID}CpCPREPLACE.R
