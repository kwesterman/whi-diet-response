#!/bin/bash


DIR=$1


## Clean genotypes
#sbatch --mem 10G --wrap="python clean_phenotypes/clean_phenotypes_bprhs.py"
#sbatch --mem 10G --wrap="python clean_phenotypes/clean_phenotypes_fhs.py"
#sbatch --mem 10G --wrap="python clean_phenotypes/clean_phenotypes_mesa.py"
#sbatch --mem 10G --wrap="python clean_phenotypes/clean_phenotypes_whi.py"

# Run GWAS
sbatch --mem 30G -t 20:00 run_gwas_int/run_gwas_bprhs.sh
sbatch --mem 60G -t 45:00 run_gwas_int/run_gwas_fhs.sh
cat mesa_races.txt | while read line; do sbatch --mem 50G -t 30:00 --wrap="run_gwas_int/run_gwas_mesa.sh $line"; done
cat whi_races.txt | while read line; do sbatch --mem 50G -t 45:00 --wrap="run_gwas_int/run_gwas_whi.sh $line"; done

# Meta-analysis
sbatch --mem 20G -t 30:00 --begin=now+45minutes --wrap="./meta_analysis.sh $DIR"

# Post-GWAS (reformatting of M-A results + QQ and Manhattan plots)
sbatch --mem 40G -t 20:00 --begin=now+65minutes --wrap="python post_gwas.py $DIR/all_meta.meta"

# Generate scores in WHI
sbatch --mem 60G -t 1:00:00 --begin=now+80minutes --wrap="./pt_model.sh $DIR"
