#!/bin/bash


module load plink2

RACE=$1
GENODIR=../data/processed/mesa
PHENOFILE=$GENODIR/mesa_${RACE}_gwas_phenos_sbp.txt
OUT_PREFIX=../data/processed/sfaR_sbpR/mesa_$RACE
PHENO=sfaR_sbpR

plink2 --pfile $GENODIR/mesa_$RACE \
	--pheno $PHENOFILE \
	--pheno-name ${PHENO}_INT \
	--update-sex $PHENOFILE \
	--glm sex a0-ref \
	--out $OUT_PREFIX
	#--covar-name \

head -1 $OUT_PREFIX.${PHENO}_INT.glm.linear > $OUT_PREFIX.res
awk '$7 == "ADD"' $OUT_PREFIX.${PHENO}_INT.glm.linear >> $OUT_PREFIX.res
