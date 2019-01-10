#!/bin/bash


module load plink2

GENODIR=../data/processed/fhs
PHENOFILE=$GENODIR/fhs_gwas_phenos.txt
OUT_PREFIX=../data/processed/f2c_tg/fhs

plink2 --pfile $GENODIR/fhs \
	--pheno $PHENOFILE \
	--pheno-name f2c_tg_INT_5 \
	--update-sex $PHENOFILE \
	--covar-name age_5 \
	--glm sex a0-ref \
	--out $OUT_PREFIX

head -1 $OUT_PREFIX.f2c_tg_INT_5.glm.linear > $OUT_PREFIX.res
awk '$7 == "ADD"' $OUT_PREFIX.f2c_tg_INT_5.glm.linear >> $OUT_PREFIX.res
