#!/bin/bash


RF=$1
TRF=$2
PREFIX=fat_pct_binary_${TRF}_whi_white
SS=${PREFIX}.res_annot
NOMINAL_PREFIX=${PREFIX}_nominalME
NOMINAL_SS=${NOMINAL_PREFIX}.res_annot

head -1 $SS > $NOMINAL_SS
awk 'NR == FNR { a[$1]; next } $1 in a' ../whi_subsets/${RF}_nominal_snps.txt $SS >> $NOMINAL_SS

plink --bfile ../whi/whi_white_DM \
	--clump $NOMINAL_SS \
	--clump-field P \
	--clump-p1 0.05 \
	--out $NOMINAL_PREFIX

awk 'NR == FNR { a[$3]; next } $1 in a' $NOMINAL_PREFIX.clumped $NOMINAL_SS > $NOMINAL_PREFIX.weights

plink --bfile ../whi/whi_DM \
	--score $NOMINAL_PREFIX.weights 1 6 9 \
	--out $NOMINAL_PREFIX
