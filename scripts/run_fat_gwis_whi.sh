#!/bin/bash


module load plink2

GENODIR=../whi
RACE=white
DIETVAR=fat_pct_binary
PHENO=$1
PHENOFILE=whi_${RACE}_gwas_phenos.txt
OUT_PREFIX=${DIETVAR}_${PHENO}_whi_${RACE}

plink2 --pfile $GENODIR/whi \
	--covar $PHENOFILE \
	--split-cat-pheno omit-last covar-01 dataset \
	--covar-name $DIETVAR age tot_cal PC1 PC2 PC3 PC4 PC5 dataset \
	--write-covar \
	--out $OUT_PREFIX

plink2 --pfile $GENODIR/whi \
	--pheno $PHENOFILE \
	--pheno-name $PHENO \
	--update-sex $PHENOFILE \
	--covar $OUT_PREFIX.cov \
	--glm interaction a0-ref \
	--parameters 1-14,22 \
	--write-covar \
	--out $OUT_PREFIX

head -1 $OUT_PREFIX.${PHENO}.glm.linear > $OUT_PREFIX.res
awk -v dv="$DIETVAR" '$7 == "ADDx"dv' $OUT_PREFIX.${PHENO}.glm.linear >> $OUT_PREFIX.res
