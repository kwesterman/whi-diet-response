#!/bin/bash


PHENO=$1
DIR=$2


## Clean phenotypes
#sbatch --mem 5G --wrap="python clean_phenotypes_statinAdj/clean_phenotypes_bprhs.py"
#sbatch --mem 5G --wrap="python clean_phenotypes_statinAdj/clean_phenotypes_fhs.py"
#sbatch --mem 5G --wrap="python clean_phenotypes_statinAdj/clean_phenotypes_mesa.py"
#sbatch --mem 5G --wrap="python clean_phenotypes_statinAdj/clean_phenotypes_whi.py"

# Run GWAS
#sbatch --mem 30G -t 15:00 run_gwas_statinAdj/run_gwas_bprhs.sh
sbatch --mem 70G -t 90:00 run_gwas_statinAdj/run_gwas_fhs.sh $PHENO $DIR
cat mesa_races.txt | while read line; do sbatch --mem 50G -t 45:00 --wrap="run_gwas_statinAdj/run_gwas_mesa.sh $line $PHENO $DIR"; done
cat whi_races.txt | while read line; do sbatch --mem 80G -t 120:00 --wrap="run_gwas_statinAdj/run_gwas_whi.sh $line $PHENO $DIR"; done

# Meta-analysis
sbatch --mem 30G -t 40:00 --begin=now+120minutes --wrap="./meta_analysis_all.sh $DIR"
sbatch --mem 30G -t 40:00 --begin=now+120minutes --wrap="./meta_analysis_white.sh $DIR"
sbatch --mem 30G -t 40:00 --begin=now+120minutes --wrap="./meta_analysis_whi.sh $DIR"

## Calculate scores in WHI
#sbatch --mem 50G -t 40:00 --begin=now+90minutes --wrap="./pt_model.sh $DIR"
