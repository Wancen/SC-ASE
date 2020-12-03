#!/bin/bash

#SBATCH -p general
#SBATCH -N 1
#SBATCH --mem=16g
#SBATCH -n 6
#SBATCH -t 2:00:00

R CMD BATCH --no-save --no-restore '--args 10 50 csv/test.csv' smurf.R log/test.Rout
