#!/bin/bash


module load plink2

RACE=$1
GENODIR=../data/processed/whi
PHENOFILE=$GENODIR/whi_${RACE}_gwas_phenos_sbp.txt
OUT_PREFIX=../data/processed/sfaR_sbpR/whi_$RACE
PHENO=sfaR_sbpR

plink2 --pfile $GENODIR/whi \
	--pheno $PHENOFILE \
	--pheno-name ${PHENO}_INT \
	--update-sex $PHENOFILE \
	--covar-name PC1 PC2 PC3 PC4 PC5 \
	--glm a0-ref \
	--out $OUT_PREFIX

head -1 $OUT_PREFIX.${PHENO}_INT.glm.linear > $OUT_PREFIX.res
awk '$7 == "ADD"' $OUT_PREFIX.${PHENO}_INT.glm.linear >> $OUT_PREFIX.res