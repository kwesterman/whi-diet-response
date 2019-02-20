#!/bin/bash


module load plink2
module load gcta


GENODIR=../data/processed/whi
DIR=../data/processed/architecture

plink2 --pfile $GENODIR/whi --keep $DIR/whi_white_cs_phenos.txt --maf 0.01 --hard-call-threshold 0 --geno 0 --make-grm-bin --out $DIR/whi_white_cs
plink2 --pfile $GENODIR/whi --keep $DIR/whi_white_delta_phenos.txt --maf 0.01 --hard-call-threshold 0 --geno 0 --make-grm-bin --out $DIR/whi_white_delta
plink2 --pfile $GENODIR/whi --keep $DIR/whi_white_trial_phenos.txt --maf 0.01 --hard-call-threshold 0 --geno 0 --make-grm-bin --out $DIR/whi_white_trial

python << EOF
import pandas as pd

phenos = pd.read_csv("$DIR/whi_white_" + vartype + "_phenos.txt", sep=" ")
phenos.filter(["FID", "IID", "logf2c_bmi_INT"]).to_csv("$DIR/whi_cs_logf2c_bmi.txt", sep=" ", index=False)
phenos.filter(["FID", "IID", "logf2c_sbp_INT"]).to_csv("$DIR/whi_cs_logf2c_sbp.txt", sep=" ", index=False)

phenos = pd.read_csv("$DIR/whi_white_" + vartype + "_phenos.txt", sep=" ")
phenos.filter(["FID", "IID", "logf2c_bmi_INT"]).to_csv("$DIR/whi_cs_logf2c_bmi.txt", sep=" ", index=False)
phenos.filter(["FID", "IID", "logf2c_sbp_INT"]).to_csv("$DIR/whi_cs_logf2c_sbp.txt", sep=" ", index=False)

phenos = pd.read_csv("$DIR/whi_white_" + vartype + "_phenos.txt", sep=" ")
phenos.filter(["FID", "IID", "logf2c_bmi_INT"]).to_csv("$DIR/whi_cs_logf2c_bmi.txt", sep=" ", index=False)
phenos.filter(["FID", "IID"] + ["PC" + str(idx) for idx in range(1,6)]).to_csv("$DIR/whi_cs_PCs.txt", sep=" ", index=False)
EOF

gcta64 --reml --grm-bin $DIR/whi_cs --pheno $DIR/whi_cs_logf2c_bmi.txt --qcovar $DIR/whi_cs_PCs.txt --out $DIR/whi_cs_bmi
