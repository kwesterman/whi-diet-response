#!/bin/bash


RF=$1
TRF=$2
PREFIX=fat_pct_binary_${TRF}_whi_white
SS=${PREFIX}.res_annot
ALL_PREFIX=${PREFIX}_all
# ALL_SS=${ALL_PREFIX}.res_annot

# head -1 $SS > $ALL_SS
# awk 'NR == FNR { a[$1]; next } $1 in a' ../whi_subsets/${RF}_all_snps.txt $SS >> $ALL_SS

plink --bfile ../whi/whi_white_DM \
	--clump $SS \
	--clump-field P \
	--clump-p1 0.05 \
	--out $ALL_PREFIX

awk 'NR == FNR { a[$3]; next } $1 in a' $ALL_PREFIX.clumped $SS > $ALL_PREFIX.weights

plink --bfile ../whi/whi_DM \
	--score $ALL_PREFIX.weights 1 6 9 \
	--out $ALL_PREFIX
