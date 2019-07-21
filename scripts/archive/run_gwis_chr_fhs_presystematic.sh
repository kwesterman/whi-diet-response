#!/bin/bash

module load python/3.6.0
module load plink

# Properly format phenotype and covariate data
python << EOF
import pandas as pd
import numpy as np
phenos = (pd.read_csv("../int/metaData_fhs.csv")
	.assign(sfa_pct_5 = lambda x: x.sfa_5 * 9 / x.tot_cal_5)
	.query('lipid_med5 == False')
	.loc[:,['subjID','subjID','ldl5','sfa_pct_5','age5','sex']]
	.dropna()
)
phenos.columns = ['FID','IID'] + list(phenos.columns)[2:]
phenos.to_csv("../int/plink_phenos_fhs.txt", sep=" ", na_rep='NA', index=False)
EOF

CHR=$SLURM_ARRAY_TASK_ID
OUT_PREFIX=../int/gwis/fhs_gwis_res_chr$CHR

echo "$CHR 1 500000000 chr$CHR" > extractChr${CHR}.txt

plink --bfile ../int/plink_imputed/fhs_genos \
--extract range extractChr${CHR}.txt \
--pheno ../int/plink_phenos_fhs.txt --pheno-name ldl5 \
--update-sex ../int/plink_phenos_fhs.txt 4 \
--covar ../int/plink_phenos_fhs.txt --covar-name sfa_pct_5 age5 \
--linear sex interaction --parameters 1-4,6 \
--ci 0.95 \
--out $OUT_PREFIX

awk '$5~/TEST|ADDxsfa_pct_5/' ${OUT_PREFIX}.assoc.linear > ${OUT_PREFIX}_int.txt

plink --bfile ../int/plink_imputed/fhs_genos \
--clump ${OUT_PREFIX}_int.txt \
--out $OUT_PREFIX

#grep -E 'SNP|rs' ../int/gwis/fhs_gwis_res.assoc.linear > ../int/gwis/fhs_gwis_res.txt
