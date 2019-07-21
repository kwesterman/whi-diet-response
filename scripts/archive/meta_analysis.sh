#!/bin/bash


module load plink

DIR=$1

# Meta-analysis of all relevant cohorts
declare -a gwas=("whi_white" "whi_black" "whi_hispanic" "mesa_white" "mesa_black" "mesa_hispanic" "mesa_asian")
gwas_files=("${gwas[@]/%/.res}")
gwas_filepaths=("${gwas_files[@]/#/$DIR/}")

plink --meta-analysis ${gwas_filepaths[@]} + qt no-map \
	--meta-analysis-snp-field ID \
	--out $DIR/all_meta

# Meta-analysis of white subgroups only
declare -a gwas_white=("fhs" "whi_white" "mesa_white")
gwas_white_files=("${gwas_white[@]/%/.res}")
gwas_white_filepaths=("${gwas_white_files[@]/#/$DIR/}")

plink --meta-analysis ${gwas_white_filepaths[@]} + qt no-map \
	--meta-analysis-snp-field ID \
	--out $DIR/white_meta	

## Meta-analysis of non-white subgroups only
#declare -a gwas_nonwhite=("whi_black" "whi_hispanic" "mesa_black" "mesa_hispanic" "mesa_asian" "bprhs")
#gwas_nonwhite_files=("${gwas_nonwhite[@]/%/.res}")
#gwas_nonwhite_filepaths=("${gwas_nonwhite_files[@]/#/$DIR/}")
#
#plink --meta-analysis ${gwas_nonwhite_filepaths[@]} + qt no-map \
#	--meta-analysis-snp-field ID \
#	--out $DIR/nonwhite_meta	
#
## Meta-analysis of WHI only
#declare -a gwas_whi=("whi_white" "whi_black" "whi_hispanic")
#gwas_whi_files=("${gwas_whi[@]/%/.res}")
#gwas_whi_filepaths=("${gwas_whi_files[@]/#/$DIR/}")
#
#plink --meta-analysis ${gwas_whi_filepaths[@]} + qt no-map \
#	--meta-analysis-snp-field ID \
#	--out $DIR/whi_meta	
#
## Meta-analysis of MESA only
#declare -a gwas_mesa=("mesa_white" "mesa_black" "mesa_hispanic" "mesa_asian")
#gwas_mesa_files=("${gwas_mesa[@]/%/.res}")
#gwas_mesa_filepaths=("${gwas_mesa_files[@]/#/$DIR/}")
#
#plink --meta-analysis ${gwas_mesa_filepaths[@]} + qt no-map \
#	--meta-analysis-snp-field ID \
#	--out $DIR/mesa_meta	
