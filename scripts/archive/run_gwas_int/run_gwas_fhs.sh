#!/bin/bash


module load plink2

GENODIR=../data/processed/fhs
PHENOFILE=$GENODIR/fhs_gwas_phenos_long.txt
DIET=$1
PHENO=$2
DIR=$3
OUT_PREFIX=$DIR/fhs

plink2 --pfile $GENODIR/fhs \
	--pheno $PHENOFILE \
	--pheno-name ${PHENO}_INT \
	--update-sex $PHENOFILE \
	--covar-name $DIET age \
	--glm sex interaction a0-ref \
	--parameters 1-3,5,7 \
	--out $OUT_PREFIX

head -1 $OUT_PREFIX.${PHENO}_INT.glm.linear > $OUT_PREFIX.main.res
awk '$7 == "ADD"' $OUT_PREFIX.${PHENO}_INT.glm.linear >> $OUT_PREFIX.main.res
head -1 $OUT_PREFIX.${PHENO}_INT.glm.linear > $OUT_PREFIX.interaction.res
awk -v diet=$DIET '$7 == "ADDx"diet' $OUT_PREFIX.${PHENO}_INT.glm.linear >> $OUT_PREFIX.interaction.res
