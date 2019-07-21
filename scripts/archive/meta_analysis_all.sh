#!/bin/bash


module load plink

DIR=$1

# Meta-analysis of all relevant cohorts
declare -a gwas=("whi_white" "whi_black" "whi_hispanic" "fhs" "mesa_white" "mesa_black" "mesa_hispanic" "mesa_asian")
gwas_files=("${gwas[@]/%/.res}")
gwas_filepaths=("${gwas_files[@]/#/$DIR/}")

plink --meta-analysis ${gwas_filepaths[@]} + qt no-map \
	--meta-analysis-snp-field ID \
	--out $DIR/all_meta

python post_gwas.py $DIR/all_meta.meta
