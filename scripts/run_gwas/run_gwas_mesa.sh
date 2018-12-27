#!/bin/bash


module load plink2

RACE=$1
GENODIR=../data/processed/mesa
PHENOFILE=$GENODIR/mesa_${RACE}_gwas_phenos.txt
OUT_PREFIX=../data/processed/sfa_ldl/mesa_$RACE

plink2 --pfile $GENODIR/mesa_$RACE \
	--pheno $PHENOFILE \
	--pheno-name sfa_ldl_product \
	--update-sex $PHENOFILE \
	--covar-name age bmi pufa_pct \
	--glm sex a0-ref \
	--out $OUT_PREFIX

head -1 $OUT_PREFIX.sfa_ldl_product.glm.linear > $OUT_PREFIX.res
awk '$7 == "ADD"' $OUT_PREFIX.sfa_ldl_product.glm.linear >> $OUT_PREFIX.res
