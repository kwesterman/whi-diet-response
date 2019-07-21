#!/bin/bash


module load plink2

GENODIR=../data/processed/fhs
PHENOFILE=$GENODIR/fhs_gwas_phenos_long.txt
PHENO=$1
DIR=$2
OUT_PREFIX=$DIR/fhs

plink2 --pfile $GENODIR/fhs \
	--pheno $PHENOFILE \
	--pheno-name ${PHENO}_INT \
	--update-sex $PHENOFILE \
	--glm sex a0-ref \
	--out $OUT_PREFIX
	#--covar-name  \

head -1 $OUT_PREFIX.${PHENO}_INT.glm.linear > $OUT_PREFIX.res
awk '$7 == "ADD"' $OUT_PREFIX.${PHENO}_INT.glm.linear >> $OUT_PREFIX.res
