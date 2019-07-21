#!/bin/bash


DIET=$1
PHENO=$2
DIR=$3

mkdir -p $DIR

# Run GWAS
sbatch --mem 60G -t 5:00:00 run_gwas_int/run_gwas_fhs.sh $DIET $PHENO $DIR
sbatch --mem 80G -t 2:00:00 run_gwas_int/run_gwas_whi.sh white $DIET $PHENO $DIR
sbatch --mem 50G -t 2:00:00 run_gwas_int/run_gwas_mesa.sh white $DIET $PHENO $DIR

# Meta-analysis
sbatch --mem 30G -t 40:00 --begin=now+5hours --wrap="./meta_analysis_white_int.sh $DIR"

## Calculate scores in WHI
#sbatch --mem 50G -t 40:00 --begin=now+90minutes --wrap="./pt_model.sh $DIR"
