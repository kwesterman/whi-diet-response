#!/bin/bash


module load plink2

GENODIR=../data/processed/fhs
PHENOFILE=$GENODIR/fhs_gwas_phenos.txt
OUT_PREFIX=../data/processed/sfa_ldl_interaction/fhs

plink2 --pfile $GENODIR/fhs \
	--pheno $PHENOFILE \
	--pheno-name ldl_5 \
	--update-sex $PHENOFILE \
	--covar-name sfa_pct_5 age_5 \
	--glm sex interaction a0-ref \
	--parameters 1-4, 7 \
	--out $OUT_PREFIX

head -1 $OUT_PREFIX.ldl_5.glm.linear > $OUT_PREFIX.res
awk '$7 == "ADDxsfa_pct_5"' $OUT_PREFIX.ldl_5.glm.linear >> $OUT_PREFIX.res
