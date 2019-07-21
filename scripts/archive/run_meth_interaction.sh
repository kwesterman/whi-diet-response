#!/bin/bash

#SBATCH -p batch,largemem
#SBATCH --mem 60G
#SBATCH -c 1
#SBATCH -t 4:00:00
#SBATCH -o meth_interaction-%A.out

module load R/3.4.3
time Rscript -e "Sys.setenv(RSTUDIO_PANDOC='/opt/shared/RStudio/0.98/bin/pandoc'); "\
"rmarkdown::render('meth_interaction.Rmd', output_dir='../doc/meth_interaction/')"
