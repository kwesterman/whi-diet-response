#!/bin/bash

module load plink

DIR=$1

CHR=$SLURM_ARRAY_TASK_ID
OUT_PREFIX=$DIR/mesa_cau_res_chr$CHR
PHENO=int_product

echo "$CHR 1 500000000 chr$CHR" > $DIR/extractChr${CHR}.txt

plink --bfile ../int/mesa/mesa_cau \
	--allow-extra-chr \
	--pheno $DIR/mesa_cau_phenos.txt --pheno-name $PHENO \
	--chr $CHR \
	--update-sex $DIR/mesa_cau_phenos.txt \
	--keep $DIR/mesa_cau_phenos.txt \
	--covar $DIR/mesa_cau_phenos.txt \
	--covar-name age bmi pufa_pct \
	--linear sex --parameters 1-5 \
	--ci 0.95 \
	--out $OUT_PREFIX

head -1 $OUT_PREFIX.assoc.linear > $OUT_PREFIX.txt
awk '$5 == "ADD"' $OUT_PREFIX.assoc.linear >> $OUT_PREFIX.txt
echo $CHR >> $DIR/mesa_cau_complete.txt
