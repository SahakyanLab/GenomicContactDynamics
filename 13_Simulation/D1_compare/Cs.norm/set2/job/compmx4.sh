#!/bin/sh
##########################################################################
## A script template for submitting batch jobs.
## Please note that anything after "#SBATCH" on a line will be treated as
## a SLURM command.
##########################################################################
#SBATCH --mem=18G
#SBATCH -n 4
#SBATCH --array=1-21
#SBATCH --mail-user=ltamon
#SBATCH --mail-type=ALL
##########################################################################
module load R-base/4.0.1
module load R-cbrg/current

R --vanilla < /t1-data/user/ltamon/DPhil/GenomicContactDynamics/21_Simulation/B1_compare/Cs.norm/set2/script/ct${SLURM_ARRAY_TASK_ID}subj4ref1.R
