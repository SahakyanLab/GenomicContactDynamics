#!/bin/sh
##########################################################################
## A script template for submitting batch jobs.
## Please note that anything after "#SBATCH" on a line will be treated as
## a SLURM command.
##########################################################################
#SBATCH --mem=2G
#SBATCH -n 1
#SBATCH --array=1-7
#SBATCH --mail-user=ltamon
#SBATCH --mail-type=ALL
##########################################################################
module load R-base/4.0.1
module load R-cbrg/current

R --vanilla < /t1-data/user/ltamon/DPhil/GenomicContactDynamics/19_Mutation_rates/A2_plotsVsCp/script/calcTmutmut${SLURM_ARRAY_TASK_ID}.R
