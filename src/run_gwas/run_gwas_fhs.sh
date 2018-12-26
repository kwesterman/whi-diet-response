#!/bin/bash


module load plink2

GENODIR=../data/processed/fhs
PHENOFILE=$GENODIR/fhs_gwas_phenos.txt
OUT_PREFIX=../data/processed/sfa_ldl/fhs

plink2 --pfile $GENODIR/fhs \
	--pheno $PHENOFILE \
	--pheno-name sfa_ldl_product_5 \
	--update-sex $PHENOFILE \
	--covar-name age_5 bmi_5 pufa_pct_5 \
	--glm sex a0-ref \
	--out $OUT_PREFIX

head -1 $OUT_PREFIX.sfa_ldl_product_5.glm.linear > $OUT_PREFIX.res
awk '$7 == "ADD"' $OUT_PREFIX.sfa_ldl_product_5.glm.linear >> $OUT_PREFIX.res
