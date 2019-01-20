#!/bin/bash


module load plink2

RACE=$1
GENODIR=../data/processed/mesa
PHENOFILE=$GENODIR/mesa_${RACE}_gwas_phenos.txt
OUT_PREFIX=../data/processed/sfa_ldl_interaction/mesa_$RACE

plink2 --pfile $GENODIR/mesa_$RACE \
	--pheno $PHENOFILE \
	--pheno-name ldl \
	--update-sex $PHENOFILE \
	--covar-name sfa_pct age \
	--glm sex interaction a0-ref \
	--parameters 1-4, 7 \
	--out $OUT_PREFIX

head -1 $OUT_PREFIX.ldl.glm.linear > $OUT_PREFIX.res
awk '$7 == "ADDxsfa_pct"' $OUT_PREFIX.ldl.glm.linear >> $OUT_PREFIX.res
