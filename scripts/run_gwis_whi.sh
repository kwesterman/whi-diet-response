#!/bin/bash


module load plink

GENODIR=../data/processed/whi
RACE=white
PHENO=$1
PHENOFILE=whi_white_gwas_phenos.txt
OUT_PREFIX=../data/processed/gen6/${PHENO}_whi_$RACE

plink --bfile $GENODIR/whi_hardcalls \
	--maf 0.01 \
	--pheno $PHENOFILE \
	--pheno-name $PHENO \
	--update-sex $PHENOFILE \
	--covar $PHENOFILE \
	--covar-name fat_binary age PC1 PC2 PC3 PC4 PC5 \
	--linear interaction \
	--parameters 1-8,15 \
	--out $OUT_PREFIX
	

head -1 $OUT_PREFIX.assoc.linear > fat_$OUT_PREFIX.res
awk '$5 == "ADDxfat_binary"' $OUT_PREFIX.assoc.linear >> fat_$OUT_PREFIX.res
