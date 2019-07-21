#!/bin/bash


module load plink2
module load gcta


GENODIR=../data/processed/whi
DIR=../data/processed/architecture

### WHITES ONLY ###
#plink2 --pfile $GENODIR/whi --keep $DIR/whi_white_cs_phenos.txt --maf 0.01 --hard-call-threshold 0 --geno 0 --make-grm-bin --out $DIR/whi_white_cs

python << EOF
import pandas as pd

phenos = pd.read_csv("$DIR/whi_white_cs_phenos.txt", sep=" ")
phenos.filter(["FID", "IID", "logf2c_bmi_INT"]).to_csv("$DIR/whi_white_cs_logf2c_bmi.txt", sep=" ", index=False)
phenos.filter(["FID", "IID", "logf2c_sbp_INT"]).to_csv("$DIR/whi_white_cs_logf2c_sbp.txt", sep=" ", index=False)
phenos.filter(["FID", "IID", "sfa_bmi_INT"]).to_csv("$DIR/whi_white_cs_sfa_bmi.txt", sep=" ", index=False)
phenos.filter(["FID", "IID", "sfa_sbp_INT"]).to_csv("$DIR/whi_white_cs_sfa_sbp.txt", sep=" ", index=False)
phenos.filter(["FID", "IID", "sodium_sbp_INT"]).to_csv("$DIR/whi_white_cs_sodium_sbp.txt", sep=" ", index=False)
phenos.filter(["FID", "IID"] + ["PC" + str(idx) for idx in range(1,6)]).to_csv("$DIR/whi_white_cs_PCs.txt", sep=" ", index=False)
EOF

#gcta64 --reml --grm-bin $DIR/whi_white_cs --pheno $DIR/whi_white_cs_logf2c_bmi.txt --qcovar $DIR/whi_white_cs_PCs.txt --out $DIR/whi_white_cs_logf2c_bmi
#gcta64 --reml --grm-bin $DIR/whi_white_cs --pheno $DIR/whi_white_cs_logf2c_sbp.txt --qcovar $DIR/whi_white_cs_PCs.txt --out $DIR/whi_white_cs_logf2c_sbp
#gcta64 --reml --grm-bin $DIR/whi_white_cs --pheno $DIR/whi_white_cs_sfa_bmi.txt --qcovar $DIR/whi_white_cs_PCs.txt --out $DIR/whi_white_cs_sfa_bmi
#gcta64 --reml --grm-bin $DIR/whi_white_cs --pheno $DIR/whi_white_cs_sfa_sbp.txt --qcovar $DIR/whi_white_cs_PCs.txt --out $DIR/whi_white_cs_sfa_sbp
gcta64 --reml --grm-bin $DIR/whi_white_cs --pheno $DIR/whi_white_cs_sodium_sbp.txt --qcovar $DIR/whi_white_cs_PCs.txt --out $DIR/whi_white_cs_sodium_sbp

### BLACKS ONLY ###
#plink2 --pfile $GENODIR/whi --keep $DIR/whi_black_cs_phenos.txt --maf 0.01 --hard-call-threshold 0 --geno 0 --make-grm-bin --out $DIR/whi_black_cs

python << EOF
import pandas as pd

phenos = pd.read_csv("$DIR/whi_black_cs_phenos.txt", sep=" ")
phenos.filter(["FID", "IID", "logf2c_bmi_INT"]).to_csv("$DIR/whi_black_cs_logf2c_bmi.txt", sep=" ", index=False)
phenos.filter(["FID", "IID", "logf2c_sbp_INT"]).to_csv("$DIR/whi_black_cs_logf2c_sbp.txt", sep=" ", index=False)
phenos.filter(["FID", "IID", "sfa_bmi_INT"]).to_csv("$DIR/whi_black_cs_sfa_bmi.txt", sep=" ", index=False)
phenos.filter(["FID", "IID", "sfa_sbp_INT"]).to_csv("$DIR/whi_black_cs_sfa_sbp.txt", sep=" ", index=False)
phenos.filter(["FID", "IID", "sodium_sbp_INT"]).to_csv("$DIR/whi_black_cs_sodium_sbp.txt", sep=" ", index=False)
phenos.filter(["FID", "IID"] + ["PC" + str(idx) for idx in range(1,6)]).to_csv("$DIR/whi_black_cs_PCs.txt", sep=" ", index=False)
EOF

#gcta64 --reml --grm-bin $DIR/whi_black_cs --pheno $DIR/whi_black_cs_logf2c_bmi.txt --qcovar $DIR/whi_black_cs_PCs.txt --out $DIR/whi_black_cs_logf2c_bmi
#gcta64 --reml --grm-bin $DIR/whi_black_cs --pheno $DIR/whi_black_cs_logf2c_sbp.txt --qcovar $DIR/whi_black_cs_PCs.txt --out $DIR/whi_black_cs_logf2c_sbp
#gcta64 --reml --grm-bin $DIR/whi_black_cs --pheno $DIR/whi_black_cs_sfa_bmi.txt --qcovar $DIR/whi_black_cs_PCs.txt --out $DIR/whi_black_cs_sfa_bmi
#gcta64 --reml --grm-bin $DIR/whi_black_cs --pheno $DIR/whi_black_cs_sfa_sbp.txt --qcovar $DIR/whi_black_cs_PCs.txt --out $DIR/whi_black_cs_sfa_sbp
gcta64 --reml --grm-bin $DIR/whi_black_cs --pheno $DIR/whi_black_cs_sodium_sbp.txt --qcovar $DIR/whi_black_cs_PCs.txt --out $DIR/whi_black_cs_sodium_sbp


#### ALL ANCESTRIES ###
#plink2 --pfile $GENODIR/whi --keep $DIR/whi_all_cs_phenos.txt --maf 0.01 --hard-call-threshold 0 --geno 0 --make-grm-bin --out $DIR/whi_all_cs
#
#python << EOF
#import pandas as pd
#
#phenos = pd.read_csv("$DIR/whi_all_cs_phenos.txt", sep=" ")
#phenos.filter(["FID", "IID", "logf2c_bmi_INT"]).to_csv("$DIR/whi_all_cs_logf2c_bmi.txt", sep=" ", index=False)
#phenos.filter(["FID", "IID", "logf2c_sbp_INT"]).to_csv("$DIR/whi_all_cs_logf2c_sbp.txt", sep=" ", index=False)
#phenos.filter(["FID", "IID", "sfa_bmi_INT"]).to_csv("$DIR/whi_all_cs_sfa_bmi.txt", sep=" ", index=False)
#phenos.filter(["FID", "IID", "sfa_sbp_INT"]).to_csv("$DIR/whi_all_cs_sfa_sbp.txt", sep=" ", index=False)
#phenos.filter(["FID", "IID"] + ["PC" + str(idx) for idx in range(1,6)]).to_csv("$DIR/whi_all_cs_PCs.txt", sep=" ", index=False)
#EOF
#
##gcta64 --reml --grm-bin $DIR/whi_all_cs --pheno $DIR/whi_all_cs_logf2c_bmi.txt --qcovar $DIR/whi_all_cs_PCs.txt --out $DIR/whi_all_cs_logf2c_bmi
##gcta64 --reml --grm-bin $DIR/whi_all_cs --pheno $DIR/whi_all_cs_logf2c_sbp.txt --qcovar $DIR/whi_all_cs_PCs.txt --out $DIR/whi_all_cs_logf2c_sbp
#gcta64 --reml --grm-bin $DIR/whi_all_cs --pheno $DIR/whi_all_cs_sfa_bmi.txt --qcovar $DIR/whi_all_cs_PCs.txt --out $DIR/whi_all_cs_sfa_bmi
#gcta64 --reml --grm-bin $DIR/whi_all_cs --pheno $DIR/whi_all_cs_sfa_sbp.txt --qcovar $DIR/whi_all_cs_PCs.txt --out $DIR/whi_all_cs_sfa_sbp
