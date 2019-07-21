#!/bin/bash


module load plink2
module load gcta

GENODIR=../data/processed/fhs
PHENOFILE=$GENODIR/fhs_gwas_phenos_long.txt
PHENO=$1
DIR=$2
OUT_PREFIX=$DIR/fhs

PHENO_COL=$(head -1 $PHENOFILE | tr ' ' '\n' | grep -xn "$PHENO" | cut -d ':' -f 1)
cut -d ' ' -f 1-3 $PHENOFILE > $GENODIR/fhs_covars.txt

gcta64 --mlma \
	--bfile $GENODIR/fhs \
	--grm $GENODIR/fhs \
	--pheno $PHENOFILE \
	--mpheno $PHENO_COL \
	--covar $GENODIR/fhs_covars.txt \
	--out $OUT_PREFIX \
	--thread-num 10

echo -e "CHR\tSNP\tBP\tALT\tREF\tFREQ\tBETA\tSE\tP" > $OUT_PREFIX.res
tail -n +2 $OUT_PREFIX.$OUT_PREFIX.mlma >> $OUT_PREFIX.res
