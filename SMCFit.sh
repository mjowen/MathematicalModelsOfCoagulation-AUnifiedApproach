#!/bin/bash
#SBATCH -p defq
#SBATCH --nodes=1
#SBATCH --ntasks-per-node 12
#SBATCH --mem=120g
#SBATCH --time=7-0


module load matlab-uon/r2021a


cd $SLURM_SUBMIT_DIR



matlab -nodisplay -nosplash < ABCSMCFitting.m
