#!/bin/bash

#SBATCH -p general
#SBATCH -N 1
#SBATCH --mem=16g
#SBATCH -n 10
#SBATCH -t 45:00

R CMD BATCH --no-save --no-restore smurf.R log/strategy2.Rout
