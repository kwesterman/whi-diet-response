#!/bin/bash


module load plink2

RACE=$1
GENODIR=../data/processed/whi
PHENOFILE=$GENODIR/whi_${RACE}_gwas_phenos_long.txt
DIET=$2
PHENO=$3
DIR=$4
OUT_PREFIX=$DIR/whi_$RACE

plink2 --pfile $GENODIR/whi \
	--pheno $PHENOFILE \
	--pheno-name ${PHENO}_INT \
	--update-sex $PHENOFILE \
	--covar-name $DIET age PC1 PC2 PC3 PC4 PC5 \
	--glm interaction a0-ref \
	--parameters 1-9 \
	--out $OUT_PREFIX

head -1 $OUT_PREFIX.${PHENO}_INT.glm.linear > $OUT_PREFIX.main.res
awk '$7 == "ADD"' $OUT_PREFIX.${PHENO}_INT.glm.linear >> $OUT_PREFIX.main.res
head -1 $OUT_PREFIX.${PHENO}_INT.glm.linear > $OUT_PREFIX.interaction.res
awk -v diet=$DIET '$7 == "ADDx"diet' $OUT_PREFIX.${PHENO}_INT.glm.linear >> $OUT_PREFIX.interaction.res
