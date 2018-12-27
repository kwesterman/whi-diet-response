#!/bin/bash


module load plink2

RACE=$1
GENODIR=../data/processed/whi
PHENOFILE=$GENODIR/whi_${RACE}_gwas_phenos.txt
OUT_PREFIX=../data/processed/sfa_ldl/whi_$RACE

plink2 --pfile $GENODIR/whi \
	--pheno $PHENOFILE \
	--pheno-name sfa_ldl_product \
	--update-sex $PHENOFILE \
	--covar-name age bmi pufa_pct PC1 PC2 PC3 PC4 PC5 \
	--glm a0-ref \
	--out $OUT_PREFIX

head -1 $OUT_PREFIX.sfa_ldl_product.glm.linear > $OUT_PREFIX.res
awk '$7 == "ADD"' $OUT_PREFIX.sfa_ldl_product.glm.linear >> $OUT_PREFIX.res
