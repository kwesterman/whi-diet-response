#!/bin/bash

module load python/3.6.0
module load plink

DIR=$1

CHR=$SLURM_ARRAY_TASK_ID
OUT_PREFIX=$DIR/whi_white_res_chr$CHR
PHENO=int_product

echo "$CHR 1 500000000 chr$CHR" > $DIR/extractChr${CHR}.txt

plink --bfile ../int/whi/whi_white_genos \
--pheno $DIR/whi_phenos.txt --pheno-name $PHENO \
--extract range $DIR/extractChr$CHR.txt \
--update-sex $DIR/whi_phenos.txt 1 \
--keep $DIR/whi_phenos.txt \
--covar $DIR/whi_phenos.txt \
--covar-name age bmi pufa PC1 PC2 PC3 PC4 PC5 \
--linear --parameters 1-9 \
--ci 0.95 \
--out $OUT_PREFIX

head -1 $OUT_PREFIX.assoc.linear > $OUT_PREFIX.txt
awk '$5 == "ADD"' $OUT_PREFIX.assoc.linear >> $OUT_PREFIX.txt
echo $CHR >> $DIR/whi_white_complete.txt

