#!/bin/bash

#SBATCH -N 1
#SBATCH --mem=16g
#SBATCH -n 10
#SBATCH -t 1:00:00

R CMD BATCH --no-save --no-restore sim2.R log/sim2_20.Rout
