#!/bin/bash

module load python/3.6.0
module load plink

DIETVAR=$1
PHENO=$2
OUTDIR=$3

python clean_metadata_whi.py $OUTDIR

CHR=$SLURM_ARRAY_TASK_ID
OUT_PREFIX=$OUTDIR/whi_res_chr$CHR

echo "$CHR 1 500000000 chr$CHR" > $OUTDIR/extractChr${CHR}.txt

if [ "$4" == "ancestry" ]; then
	COVAR_NAME="PC1 PC2 PC3 PC4 PC5 age ${DIETVAR}"
	PARAMS="1-8,15"
else
	COVAR_NAME="${DIETVAR} age"
	PARAMS="1-3,5"
fi

if [ "$4" == "whites" ]; then
	KEEP_GROUP=$OUTDIR/whi_white_ids.txt
else
	KEEP_GROUP=$OUTDIR/whi_nonDM_ids.txt
fi

plink --bfile ../int/plink_imputed/updated/whi_genos \
--pheno $OUTDIR/meta_whi.txt --pheno-name ${PHENO} \
--extract range $OUTDIR/extractChr${CHR}.txt \
--update-sex $OUTDIR/meta_whi.txt 1 \
--keep $OUTDIR/train_ids_whi.txt \
--covar $OUTDIR/meta_whi.txt --covar-name $COVAR_NAME pufa \
--linear interaction --parameters 1-9,17 \
--ci 0.95 \
--out $OUT_PREFIX

awk -v re="TEST|ADDx${DIETVAR}" '$5~re' ${OUT_PREFIX}.assoc.linear > ${OUT_PREFIX}_int.txt
echo $CHR >> $OUTDIR/whi_complete.txt

#plink --bfile ../int/plink_imputed/updated/whi_genos \
#--clump ${OUT_PREFIX}_int.txt \
#--out $OUT_PREFIX
