#!/bin/bash

module load python/3.6.0
module load plink
module load plink2

DIR=$1

CHR=$SLURM_ARRAY_TASK_ID
OUT_PREFIX=$DIR/bprhs_res_chr$CHR
PHENO=int_product

echo "$CHR 1 500000000 chr$CHR" > $DIR/extractChr${CHR}.txt

plink --bfile ../int/bprhs/bprhs_genos \
--pheno $DIR/bprhs_phenos.txt --pheno-name $PHENO \
--chr $CHR \
--update-sex $DIR/bprhs_phenos.txt \
--keep $DIR/bprhs_phenos.txt \
--covar $DIR/bprhs_phenos.txt \
--covar-name age bmi pufa PC1 \
--linear sex --parameters 1-6 \
--ci 0.95 \
--out $OUT_PREFIX

head -1 $OUT_PREFIX.assoc.linear > $OUT_PREFIX.txt
awk '$5 == "ADD"' $OUT_PREFIX.assoc.linear >> $OUT_PREFIX.txt
echo $CHR >> $DIR/bprhs_complete.txt
