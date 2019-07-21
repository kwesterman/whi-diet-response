#!/bin/bash


module load plink

DIR=$1

# Meta-analysis of WHI only
declare -a gwas_whi=("whi_white" "whi_black" "whi_hispanic")
gwas_whi_files=("${gwas_whi[@]/%/.res}")
gwas_whi_filepaths=("${gwas_whi_files[@]/#/$DIR/}")

plink --meta-analysis ${gwas_whi_filepaths[@]} + qt no-map \
	--meta-analysis-snp-field ID \
	--out $DIR/whi_meta	

python post_gwas.py $DIR/whi_meta.meta
