#!/bin/bash


RF=$1
TRF=$2
PREFIX=fat_pct_binary_${TRF}_whi_white
SS=${PREFIX}.res_annot
GARSKE_PREFIX=${PREFIX}_garske
GARSKE_SS=${GARSKE_PREFIX}.res_annot

head -1 $SS > $GARSKE_SS
awk 'NR == FNR { a[$1]; next } $1 in a' bmi_prioritized_variants.txt $SS >> $GARSKE_SS

plink --bfile ../whi/whi_white_DM \
	--clump $GARSKE_SS \
	--clump-field P \
	--clump-p1 0.05 \
	--out $GARSKE_PREFIX

awk 'NR == FNR { a[$3]; next } $1 in a' $GARSKE_PREFIX.clumped $GARSKE_SS > $GARSKE_PREFIX.weights

plink --bfile ../whi/whi_DM \
	--score $GARSKE_PREFIX.weights 1 6 9 \
	--out $GARSKE_PREFIX
