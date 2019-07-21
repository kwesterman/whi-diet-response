#!/bin/bash

#SBATCH -p batch,largemem
#SBATCH --mem 25G
#SBATCH -c 8
#SBATCH -t 3:00:00
#SBATCH -o sfa_ldl_models_functional-%A.out

module load R/3.4.3
time Rscript -e "Sys.setenv(RSTUDIO_PANDOC='/opt/shared/RStudio/0.98/bin/pandoc'); rmarkdown::render('sfa_ldl_models_functional.Rmd', output_dir='../output/')"
