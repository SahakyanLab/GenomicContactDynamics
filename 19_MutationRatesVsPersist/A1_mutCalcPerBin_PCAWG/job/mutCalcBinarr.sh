#!/bin/sh
##########################################################################
## A script template for submitting batch jobs.
## Please note that anything after "#SBATCH" on a line will be treated as
## a SLURM command.
##########################################################################
##14G 3C for All(1); 9G 1C for C>T(4); 6G 1C for C>A(2,3,5,6,7)
#SBATCH --mem=14G
#SBATCH -n 3
#SBATCH --array=1
#SBATCH --mail-user=ltamon
#SBATCH --mail-type=ALL
##########################################################################
module load R-base/4.0.1
module load R-cbrg/current

R --vanilla < /t1-data/user/ltamon/DPhil/GenomicContactDynamics/19_Mutation_rates/A1_mutCalcPerBin/script/mutCalcPerBin${SLURM_ARRAY_TASK_ID}.R
