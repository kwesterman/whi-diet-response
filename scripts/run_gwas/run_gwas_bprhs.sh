#!/bin/bash


module load plink2

GENODIR=../data/processed/bprhs
PHENOFILE=$GENODIR/bprhs_gwas_phenos.txt
OUT_PREFIX=../data/processed/sfa_ldl/bprhs

plink2 --pfile $GENODIR/bprhs \
	--pheno $PHENOFILE \
	--pheno-name sfa_ldl_product \
	--update-sex $PHENOFILE \
	--covar-name age bmi pufa_pct PC1 \
	--glm sex a0-ref \
	--out $OUT_PREFIX

head -1 $OUT_PREFIX.sfa_ldl_product.glm.linear > $OUT_PREFIX.res
awk '$7 == "ADD"' $OUT_PREFIX.sfa_ldl_product.glm.linear >> $OUT_PREFIX.res
