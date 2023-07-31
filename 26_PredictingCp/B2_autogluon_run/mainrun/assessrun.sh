#!/bin/sh
##########################################################################
## A script template for submitting batch jobs.
## Please note that anything after "#SBATCH" on a line will be treated as
## a SLURM command.
##########################################################################
#SBATCH -p batch
#SBATCH --mem=10G
#SBATCH -n 1
#SBATCH --mail-user=ltamon
#SBATCH --mail-type=ALL
##########################################################################
module load python-base/3.8.3

python /project/sahakyanlab/ltamon/SahakyanLab/GenomicContactDynamics/26_PredictingCp/B2_autogluon_run/assessrun/assessrun.py
