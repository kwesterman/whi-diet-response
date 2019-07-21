#!/bin/bash

module load python/3.6.0
module load plink

DIETVAR=$1
PHENO=$2
OUTDIR=$3
EXAM=$4
DVNAME=${DIETVAR}_${EXAM}

CHR=$SLURM_ARRAY_TASK_ID
OUT_PREFIX=$OUTDIR/fhs_res_chr$CHR

echo "$CHR 1 500000000 chr$CHR" > $OUTDIR/extractChr${CHR}.txt  # In effect: take all SNPs from given chromosome

plink --bfile ../int/plink_imputed/fhs_genos \
--pheno $OUTDIR/meta_fhs.txt --pheno-name ${PHENO}_${EXAM} \
--extract range $OUTDIR/extractChr${CHR}.txt \
--update-sex $OUTDIR/meta_fhs.txt 1 \
--keep $OUTDIR/train_ids_fhs.txt \
--covar $OUTDIR/meta_fhs.txt --covar-name $DVNAME pufa_${EXAM} age_${EXAM} \
--linear sex interaction --parameters 1-4,7-8 \
--ci 0.95 \
--out $OUT_PREFIX

awk -v re="TEST|ADDx$DVNAME" '$5~re' ${OUT_PREFIX}.assoc.linear > ${OUT_PREFIX}_int.txt
echo $CHR >> $OUTDIR/fhs_complete.txt

#plink --bfile ../int/plink_imputed/fhs_genos \
#--clump ${OUT_PREFIX}_int.txt \
#--out $OUT_PREFIX
