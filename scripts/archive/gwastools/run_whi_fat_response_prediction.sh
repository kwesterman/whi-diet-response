#!/bin/bash

#SBATCH -p batch,largemem
#SBATCH --mem 160G
#SBATCH -c 4
#SBATCH -t 6:00:00
#SBATCH -o whi_fat_response_prediction-%A.out

module load R/3.6
time Rscript -e "Sys.setenv(RSTUDIO_PANDOC='/opt/shared/RStudio/0.98/bin/pandoc'); rmarkdown::render('whi_fat_response_prediction.Rmd')"
