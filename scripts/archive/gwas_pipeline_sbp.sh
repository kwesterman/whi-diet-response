#!/bin/bash


DIR=$1


# Clean phenotypes
#sbatch --mem 5G --wrap="python clean_phenotypes/clean_phenotypes_bprhs.py"
#sbatch --mem 5G --wrap="python clean_phenotypes/clean_phenotypes_fhs.py"
#sbatch --mem 5G --wrap="python clean_phenotypes/clean_phenotypes_mesa.py"
#sbatch --mem 5G --wrap="python clean_phenotypes/clean_phenotypes_whi.py"

# Run GWAS
sbatch --mem 30G -t 15:00 run_gwas_sbp/run_gwas_bprhs.sh
sbatch --mem 60G -t 30:00 run_gwas_sbp/run_gwas_fhs.sh
cat mesa_races.txt | while read line; do sbatch --mem 50G -t 20:00 --wrap="run_gwas_sbp/run_gwas_mesa.sh $line"; done
cat whi_races.txt | while read line; do sbatch --mem 70G -t 40:00 --wrap="run_gwas_sbp/run_gwas_whi.sh $line"; done

# Meta-analysis
sbatch --mem 30G -t 30:00 --begin=now+40minutes --wrap="./meta_analysis.sh $DIR"

# Post-processing + plots
sbatch --mem 50G --begin=now+80minutes --wrap="python post_gwas.py $DIR/all_meta.meta"
sbatch --mem 50G --begin=now+80minutes --wrap="python post_gwas.py $DIR/white_meta.meta"

# Calculate scores in WHI
sbatch --mem 50G -t 30:00 --begin=now+95minutes --wrap="./pt_model.sh $DIR"
