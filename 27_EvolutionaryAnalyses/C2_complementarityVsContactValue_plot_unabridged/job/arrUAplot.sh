#!/bin/sh
##########################################################################
## A script template for submitting batch jobs.
## Please note that anything after "#SBATCH" on a line will be treated as
## a SLURM command.
##########################################################################
#SBATCH --mem=100G
#SBATCH -n 1
#SBATCH --array=1-5
#SBATCH --mail-user=ltamon
#SBATCH --mail-type=ALL
##########################################################################
module load R-base/4.2.0
module load R-cbrg/current

R --vanilla < /project/sahakyanlab/ltamon/SahakyanLab/GenomicContactDynamics/27_EvolutionaryAnalyses/C2_complementarityVsContactValue_plot_unabridged/script/UAplot${SLURM_ARRAY_TASK_ID}.R
