#!/bin/bash


PHENO=$1
DIR=$2

mkdir -p $DIR

# Run GWAS
sbatch --mem 60G -t 50:00 run_gwas/run_gwas_fhs.sh $PHENO $DIR
sbatch --mem 80G -t 60:00 run_gwas/run_gwas_whi.sh white $PHENO $DIR
sbatch --mem 50G -t 30:00 run_gwas/run_gwas_mesa.sh white $PHENO $DIR

# Meta-analysis
sbatch --mem 30G -t 40:00 --begin=now+60minutes --wrap="./meta_analysis_white.sh $DIR"

## Calculate scores in WHI
#sbatch --mem 50G -t 40:00 --begin=now+90minutes --wrap="./pt_model.sh $DIR"
