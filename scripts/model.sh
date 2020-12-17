#!/bin/bash

#SBATCH -N 1
#SBATCH --job-name=deng
#SBATCH --mem=16g
#SBATCH -n 11
#SBATCH -t 5:00:00
#SBATCH --mail-user=wancen@live.unc.edu
#SBATCH --mail-type=ALL

Rscript -e "require( 'rmarkdown' ); render('/proj/milovelab/mu/SC-ASE/scripts/Deng_s.data.Rmd', 'html_document')" --no-save --no-restore log/test.Rout
