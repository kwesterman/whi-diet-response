#!/bin/bash


PHENO=$1
DIR=$2


# Clean phenotypes
#sbatch --mem 5G --wrap="python clean_phenotypes_long/clean_phenotypes_whi.py"

# Run GWAS
cat whi_races.txt | while read line; do sbatch --mem 70G -t 40:00 --wrap="run_gwas_long/run_gwas_whi.sh $line $PHENO $DIR"; done

# Meta-analysis
sbatch --mem 30G -t 20:00 --begin=now+40minutes --wrap="./meta_analysis_long.sh $DIR"

# Post-processing + plots
sbatch --mem 50G --begin=now+60minutes --wrap="python post_gwas.py $DIR/whi_meta.meta"

# Calculate scores in WHI
sbatch --mem 50G -t 30:00 --begin=now+75minutes --wrap="./pt_model_long.sh $DIR"
