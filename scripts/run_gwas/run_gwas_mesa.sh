#!/bin/bash


module load plink2

RACE=$1
GENODIR=../data/processed/mesa
PHENOFILE=$GENODIR/mesa_${RACE}_gwas_phenos.txt
OUT_PREFIX=../data/processed/f2c_tg/mesa_$RACE

plink2 --pfile $GENODIR/mesa_$RACE \
	--pheno $PHENOFILE \
	--pheno-name f2c_tg_INT \
	--update-sex $PHENOFILE \
	--covar-name age \
	--glm sex a0-ref \
	--out $OUT_PREFIX

head -1 $OUT_PREFIX.f2c_tg_INT.glm.linear > $OUT_PREFIX.res
awk '$7 == "ADD"' $OUT_PREFIX.f2c_tg_INT.glm.linear >> $OUT_PREFIX.res
