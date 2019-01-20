#!/bin/bash


module load plink2

RACE=$1
GENODIR=../data/processed/whi
PHENOFILE=$GENODIR/whi_${RACE}_gwas_phenos.txt
OUT_PREFIX=../data/processed/sfa_ldl_interaction/whi_$RACE

plink2 --pfile $GENODIR/whi \
	--pheno $PHENOFILE \
	--pheno-name ldl \
	--update-sex $PHENOFILE \
	--covar-name sfa_pct age PC1 PC2 PC3 PC4 PC5 \
	--glm interaction a0-ref \
	--parameters 1-9 \
	--out $OUT_PREFIX

head -1 $OUT_PREFIX.ldl.glm.linear > $OUT_PREFIX.res
awk '$7 == "ADDxsfa_pct"' $OUT_PREFIX.ldl.glm.linear >> $OUT_PREFIX.res
