#!/bin/bash


PREFIX=$1

plink --bfile ../data/processed/whi/whi_white_DM \
	--clump $PREFIX.res_annot \
	--clump-snp-field SNP \
	--clump-field P \
	--clump-p1 0.0001 \
	--out $PREFIX

awk 'NR == FNR { a[$3]; next } $2 in a' $PREFIX.clumped $PREFIX.res_annot > $PREFIX.weights

plink --bfile ../data/processed/whi/whi_DM \
	--score $PREFIX.weights 2 4 7 \
	--out $PREFIX
