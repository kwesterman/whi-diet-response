#!/bin/bash


module load plink

DIR=$1

# Meta-analysis of white subgroups only
declare -a gwas_white=("fhs" "whi_white" "mesa_white")
gwas_white_files=("${gwas_white[@]/%/.res}")
gwas_white_filepaths=("${gwas_white_files[@]/#/$DIR/}")

plink --meta-analysis ${gwas_white_filepaths[@]} + qt no-map \
	--meta-analysis-snp-field ID \
	--out $DIR/white_meta	

python post_gwas.py $DIR/white_meta.meta
