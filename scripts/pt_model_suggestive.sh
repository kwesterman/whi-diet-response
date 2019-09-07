#!/bin/bash


RF=$1
TRF=$2
PREFIX=fat_pct_binary_${TRF}_whi_white
SS=${PREFIX}.res_annot
SUGGESTIVE_PREFIX=${PREFIX}_suggestiveME
SUGGESTIVE_SS=${SUGGESTIVE_PREFIX}.res_annot

head -1 $SS > $SUGGESTIVE_SS
awk 'NR == FNR { a[$1]; next } $1 in a' ../whi_subsets/${RF}_suggestive_snps.txt $SS >> $SUGGESTIVE_SS

plink --bfile ../whi/whi_white_DM \
	--clump $SUGGESTIVE_SS \
	--clump-field P \
	--clump-p1 0.05 \
	--out $SUGGESTIVE_PREFIX

awk 'NR == FNR { a[$3]; next } $1 in a' $SUGGESTIVE_PREFIX.clumped $SUGGESTIVE_SS > $SUGGESTIVE_PREFIX.weights

plink --bfile ../whi/whi_DM \
	--score $SUGGESTIVE_PREFIX.weights 1 6 9 \
	--out $SUGGESTIVE_PREFIX
