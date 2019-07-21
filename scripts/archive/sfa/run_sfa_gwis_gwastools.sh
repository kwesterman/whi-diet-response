#!/bin/bash

#SBATCH -p batch,largemem
#SBATCH --mem 40G
#SBATCH -c 8
#SBATCH -t 8:00:00
#SBATCH -o sfa_gwis_gwastools-%A.out

module load R/3.4.3
time Rscript -e "Sys.setenv(RSTUDIO_PANDOC='/opt/shared/RStudio/0.98/bin/pandoc'); rmarkdown::render('sfa_gwis_gwastools.Rmd', output_dir='../output/')"
