#!/bin/bash


module load plink2

GENODIR=../data/processed/fhs
PHENOFILE=$GENODIR/fhs_gwas_phenos_sbp.txt
OUT_PREFIX=../data/processed/sfaR_sbpR/fhs
PHENO=sfaR_sbpR

plink2 --pfile $GENODIR/fhs \
	--pheno $PHENOFILE \
	--pheno-name ${PHENO}_5_INT \
	--update-sex $PHENOFILE \
	--glm sex a0-ref \
	--out $OUT_PREFIX
	#--covar-name \

head -1 $OUT_PREFIX.${PHENO}_5_INT.glm.linear > $OUT_PREFIX.res
awk '$7 == "ADD"' $OUT_PREFIX.${PHENO}_5_INT.glm.linear >> $OUT_PREFIX.res
