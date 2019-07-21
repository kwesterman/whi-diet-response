#!/bin/bash


DIR=$1


## Clean phenotypes
#sbatch --mem 5G --wrap="python clean_phenotypes/clean_phenotypes_bprhs.py"
#sbatch --mem 5G --wrap="python clean_phenotypes/clean_phenotypes_fhs.py"
#sbatch --mem 5G --wrap="python clean_phenotypes/clean_phenotypes_mesa.py"
#sbatch --mem 5G --wrap="python clean_phenotypes/clean_phenotypes_whi.py"

# Run GWAS
#sbatch --mem 30G -t 15:00 run_gwas/run_gwas_bprhs.sh
sbatch --mem 60G -t 30:00 run_gwas/run_gwas_fhs.sh
cat mesa_races.txt | while read line; do sbatch --mem 50G -t 20:00 --wrap="run_gwas/run_gwas_mesa.sh $line"; done
cat whi_races.txt | while read line; do sbatch --mem 70G -t 40:00 --wrap="run_gwas/run_gwas_whi.sh $line"; done

# Meta-analysis
sbatch --mem 30G -t 40:00 --begin=now+40minutes --wrap="./meta_analysis_all.sh $DIR"
sbatch --mem 30G -t 40:00 --begin=now+40minutes --wrap="./meta_analysis_white.sh $DIR"
sbatch --mem 30G -t 40:00 --begin=now+40minutes --wrap="./meta_analysis_whi.sh $DIR"

## Calculate scores in WHI
#sbatch --mem 50G -t 40:00 --begin=now+90minutes --wrap="./pt_model.sh $DIR"
