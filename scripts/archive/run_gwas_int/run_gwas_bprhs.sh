#!/bin/bash


module load plink2

GENODIR=../data/processed/bprhs
PHENOFILE=$GENODIR/bprhs_gwas_phenos.txt
OUT_PREFIX=../data/processed/sfa_ldl_interaction/bprhs

plink2 --pfile $GENODIR/bprhs \
	--pheno $PHENOFILE \
	--pheno-name ldl \
	--update-sex $PHENOFILE \
	--covar-name sfa_pct age PC1 \
	--glm sex interaction a0-ref \
	--parameters 1-5, 9 \
	--out $OUT_PREFIX

head -1 $OUT_PREFIX.ldl.glm.linear > $OUT_PREFIX.res
awk '$7 == "ADDxsfa_pct"' $OUT_PREFIX.ldl.glm.linear >> $OUT_PREFIX.res
