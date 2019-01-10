#!/bin/bash


module load plink2

GENODIR=../data/processed/bprhs
PHENOFILE=$GENODIR/bprhs_gwas_phenos.txt
OUT_PREFIX=../data/processed/f2c_tg/bprhs

plink2 --pfile $GENODIR/bprhs \
	--pheno $PHENOFILE \
	--pheno-name f2c_tg_INT \
	--update-sex $PHENOFILE \
	--covar-name age PC1 \
	--glm sex a0-ref \
	--out $OUT_PREFIX

head -1 $OUT_PREFIX.f2c_tg_INT.glm.linear > $OUT_PREFIX.res
awk '$7 == "ADD"' $OUT_PREFIX.f2c_tg_INT.glm.linear >> $OUT_PREFIX.res
