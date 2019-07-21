#!/bin/bash

module load python/3.6.0
module load plink

# Properly format phenotype and covariate data
python << EOF
import pandas as pd
import numpy as np
phenos = (pd.read_csv("../int/metaData_whi.csv")
	.assign(sfa_pct = lambda x: x.sfa * 9 / x.tot_cal)
	.query('F60VTYP == 1 and lipid_med == False')
	.loc[:,['subjID','subjID','ldl','sfa_pct','age','sex']]
	.dropna()
)
phenos.columns = ['FID','IID'] + list(phenos.columns)[2:]
phenos.to_csv("../int/plink_phenos_whi.txt", sep=" ", na_rep='NA', index=False)
EOF

CHR=$SLURM_ARRAY_TASK_ID
OUT_PREFIX=../int/gwis/whi_gwis_res_chr$CHR

echo "$CHR 1 500000000 chr$CHR" > extractChr${CHR}.txt

plink --bfile ../int/plink_imputed/updated/whi_genos \
--extract range extractChr${CHR}.txt \
--pheno ../int/plink_phenos_whi.txt --pheno-name ldl \
--update-sex ../int/plink_phenos_whi.txt 4 \
--covar ../int/plink_phenos_whi.txt --covar-name sfa_pct age \
--linear interaction --parameters 1-4 \
--ci 0.95 \
--out $OUT_PREFIX

awk '$5~/TEST|ADDxsfa_pct/' ${OUT_PREFIX}.assoc.linear > ${OUT_PREFIX}_int.txt

plink --bfile ../int/plink_imputed/whi_genos \
--clump ${OUT_PREFIX}_int.txt \
--out $OUT_PREFIX
