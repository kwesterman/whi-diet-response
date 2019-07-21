#!/bin/bash


module load plink2
module load gcta


GENODIR=../data/processed/whi
DIR=../data/processed/architecture

#### WHITES ONLY ###
#plink2 --pfile $GENODIR/whi --keep $DIR/whi_white_trial_phenos.txt --maf 0.01 --hard-call-threshold 0 --geno 0 --make-grm-bin --out $DIR/whi_white_trial
#
#python << EOF
#import pandas as pd
#
#phenos = pd.read_csv("$DIR/whi_white_trial_phenos.txt", sep=" ")
#phenos.filter(["FID", "IID", "delta_bmi_INT"]).to_csv("$DIR/whi_white_trial_delta_bmi.txt", sep=" ", index=False)
#phenos.filter(["FID", "IID", "delta_sbp_INT"]).to_csv("$DIR/whi_white_trial_delta_sbp.txt", sep=" ", index=False)
#phenos.filter(["FID", "IID"] + ["PC" + str(idx) for idx in range(1,6)]).to_csv("$DIR/whi_white_trial_PCs.txt", sep=" ", index=False)
#EOF
#
#gcta64 --reml --grm-bin $DIR/whi_white_trial --pheno $DIR/whi_white_trial_delta_bmi.txt --qcovar $DIR/whi_white_trial_PCs.txt --out $DIR/whi_white_trial_bmi
#gcta64 --reml --grm-bin $DIR/whi_white_trial --pheno $DIR/whi_white_trial_delta_sbp.txt --qcovar $DIR/whi_white_trial_PCs.txt --out $DIR/whi_white_trial_sbp


### BLACKS ONLY ###
plink2 --pfile $GENODIR/whi --keep $DIR/whi_black_trial_phenos.txt --maf 0.01 --hard-call-threshold 0 --geno 0 --make-grm-bin --out $DIR/whi_black_trial

python << EOF
import pandas as pd

phenos = pd.read_csv("$DIR/whi_black_trial_phenos.txt", sep=" ")
phenos.filter(["FID", "IID", "delta_bmi_INT"]).to_csv("$DIR/whi_black_trial_delta_bmi.txt", sep=" ", index=False)
phenos.filter(["FID", "IID", "delta_sbp_INT"]).to_csv("$DIR/whi_black_trial_delta_sbp.txt", sep=" ", index=False)
phenos.filter(["FID", "IID"] + ["PC" + str(idx) for idx in range(1,6)]).to_csv("$DIR/whi_black_trial_PCs.txt", sep=" ", index=False)
EOF

gcta64 --reml --grm-bin $DIR/whi_black_trial --pheno $DIR/whi_black_trial_delta_bmi.txt --qcovar $DIR/whi_black_trial_PCs.txt --out $DIR/whi_black_trial_bmi
gcta64 --reml --grm-bin $DIR/whi_black_trial --pheno $DIR/whi_black_trial_delta_sbp.txt --qcovar $DIR/whi_black_trial_PCs.txt --out $DIR/whi_black_trial_sbp

#### ALL ETHNICITIES ###
#plink2 --pfile $GENODIR/whi --keep $DIR/whi_all_trial_phenos.txt --maf 0.01 --hard-call-threshold 0 --geno 0 --make-grm-bin --out $DIR/whi_all_trial
#
#python << EOF
#import pandas as pd
#
#phenos = pd.read_csv("$DIR/whi_all_trial_phenos.txt", sep=" ")
#phenos.filter(["FID", "IID", "delta_bmi_INT"]).to_csv("$DIR/whi_all_trial_delta_bmi.txt", sep=" ", index=False)
#phenos.filter(["FID", "IID", "delta_sbp_INT"]).to_csv("$DIR/whi_all_trial_delta_sbp.txt", sep=" ", index=False)
#phenos.filter(["FID", "IID"] + ["PC" + str(idx) for idx in range(1,6)]).to_csv("$DIR/whi_all_trial_PCs.txt", sep=" ", index=False)
#EOF
#
#gcta64 --reml --grm-bin $DIR/whi_all_trial --pheno $DIR/whi_all_trial_delta_bmi.txt --qcovar $DIR/whi_all_trial_PCs.txt --out $DIR/whi_all_trial_bmi
#gcta64 --reml --grm-bin $DIR/whi_all_trial --pheno $DIR/whi_all_trial_delta_sbp.txt --qcovar $DIR/whi_all_trial_PCs.txt --out $DIR/whi_all_trial_sbp
#
