#!/bin/bash

#SBATCH -p general
#SBATCH -N 1
#SBATCH --mem=16g
#SBATCH -n 10
#SBATCH -t 3:00:00

Rscript -e "require( 'rmarkdown' ); render('/proj/milovelab/mu/SC-ASE/scripts/Deng_s.data.Rmd', 'html_document')" --no-save --no-restore log/test.Rout
