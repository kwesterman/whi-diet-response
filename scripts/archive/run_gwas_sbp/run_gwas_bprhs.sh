#!/bin/bash


module load plink2

GENODIR=../data/processed/bprhs
PHENOFILE=$GENODIR/bprhs_gwas_phenos_sbp.txt
OUT_PREFIX=../data/processed/sfaR_sbpR/bprhs
PHENO=sfaR_sbpR

plink2 --pfile $GENODIR/bprhs \
	--pheno $PHENOFILE \
	--pheno-name ${PHENO}_INT \
	--update-sex $PHENOFILE \
	--covar-name PC1 \
	--glm sex a0-ref \
	--out $OUT_PREFIX

head -1 $OUT_PREFIX.${PHENO}_INT.glm.linear > $OUT_PREFIX.res
awk '$7 == "ADD"' $OUT_PREFIX.${PHENO}_INT.glm.linear >> $OUT_PREFIX.res
