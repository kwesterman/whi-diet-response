#!/bin/bash


module load plink2

SUMSTATS=../../raw/literature/ldl_gwas_ma_res_glgc.txt
NOMINAL_PREFIX=ldl_ME_nominalME

echo "SNP A1 BETA P" > $NOMINAL_PREFIX.sumstats
awk '$9 < 0.05 {print $3,toupper($4),$6,$9}' $SUMSTATS >> $NOMINAL_PREFIX.sumstats

plink --bfile ../whi/whi_white_DM \
	--clump $NOMINAL_PREFIX.sumstats \
	--clump-field P \
	--clump-p1 0.05 \
	--out $NOMINAL_PREFIX

awk 'NR == FNR { a[$3]; next } $1 in a' $NOMINAL_PREFIX.clumped $NOMINAL_PREFIX.sumstats > $NOMINAL_PREFIX.weights

plink2 --pfile ../whi/whi_DM \
	--score $NOMINAL_PREFIX.weights 1 2 3 \
	--out dosage_scores/$NOMINAL_PREFIX
