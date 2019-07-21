#!/bin/bash


module load plink2

GENODIR=../data/processed/bprhs
PHENOFILE=$GENODIR/bprhs_gwas_phenos_long.txt
OUT_PREFIX=../data/processed/ldl/bprhs

plink2 --pfile $GENODIR/bprhs \
	--pheno $PHENOFILE \
	--pheno-name ldl \
	--update-sex $PHENOFILE \
	--covar-name PC1 \
	--glm sex a0-ref \
	--out $OUT_PREFIX

head -1 $OUT_PREFIX.ldl.glm.linear > $OUT_PREFIX.res
awk '$7 == "ADD"' $OUT_PREFIX.ldl.glm.linear >> $OUT_PREFIX.res
