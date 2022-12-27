#!/bin/sh
##########################################################################
## A script template for submitting batch jobs.
## Please note that anything after "#SBATCH" on a line will be treated as
## SLURM command.
##########################################################################
#SBATCH --mem=60G
#SBATCH -n 1
#SBATCH --mail-user=ltamon
#SBATCH --mail-type=ALL
#########################################################################
## JOB DETAILS * JOB DETAILS * JOB DETAILS * JOB DETAILS * JOB DETAILS ##
#########################################################################
module load R-base/4.2.0
module load R-cbrg/current

Rscript --vanilla /project/sahakyanlab/ltamon/SahakyanLab/GenomicContactDynamics/27_EvolutionaryAnalyses/C2_complementarityVsContactValue_plot_unabridged.R
