#!/bin/bash

module load python/3.6.0
module load plink

DIR=$1

CHR=$SLURM_ARRAY_TASK_ID
OUT_PREFIX=$DIR/fhs_res_chr$CHR
PHENO=int_product

echo "$CHR 1 500000000 chr$CHR" > $DIR/extractChr${CHR}.txt 

plink --bfile ../int/fhs/fhs \
--pheno $DIR/fhs_phenos.txt --pheno-name ${PHENO} \
--extract range $DIR/extractChr$CHR.txt \
--update-sex $DIR/fhs_phenos.txt 1 \
--keep ../int/fhs/unrelated_ids.txt \
--covar $DIR/fhs_phenos.txt \
--covar-name age_5 bmi_5 pufa_5 \
--linear sex --parameters 1-5 \
--ci 0.95 \
--out $OUT_PREFIX

head -1 $OUT_PREFIX.assoc.linear > $OUT_PREFIX.txt
awk '$5 == "ADD"' $OUT_PREFIX.assoc.linear >> ${OUT_PREFIX}.txt
echo $CHR >> $DIR/fhs_complete.txt
