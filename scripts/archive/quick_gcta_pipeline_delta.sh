#!/bin/bash


module load plink2
module load gcta


GENODIR=../data/processed/whi
DIR=../data/processed/architecture

### WHITES ONLY ###
#plink2 --pfile $GENODIR/whi --keep $DIR/whi_white_delta_phenos.txt --maf 0.01 --hard-call-threshold 0 --geno 0 --make-grm-bin --out $DIR/whi_white_delta

python << EOF
import pandas as pd

phenos = pd.read_csv("$DIR/whi_white_delta_phenos.txt", sep=" ")
phenos.filter(["FID", "IID", "delta_logf2c_delta_bmi_INT"]).to_csv("$DIR/whi_white_delta_delta_logf2c_delta_bmi.txt", sep=" ", index=False)
phenos.filter(["FID", "IID", "delta_logf2c_delta_sbp_INT"]).to_csv("$DIR/whi_white_delta_delta_logf2c_delta_sbp.txt", sep=" ", index=False)
phenos.filter(["FID", "IID", "delta_sfa_delta_bmi_INT"]).to_csv("$DIR/whi_white_delta_delta_sfa_delta_bmi.txt", sep=" ", index=False)
phenos.filter(["FID", "IID", "delta_sfa_delta_sbp_INT"]).to_csv("$DIR/whi_white_delta_delta_sfa_delta_sbp.txt", sep=" ", index=False)
phenos.filter(["FID", "IID", "delta_sodium_delta_sbp_INT"]).to_csv("$DIR/whi_white_delta_delta_sodium_delta_sbp.txt", sep=" ", index=False)
phenos.filter(["FID", "IID"] + ["PC" + str(idx) for idx in range(1,6)]).to_csv("$DIR/whi_white_delta_PCs.txt", sep=" ", index=False)
EOF

#gcta64 --reml --grm-bin $DIR/whi_white_delta --pheno $DIR/whi_white_delta_delta_logf2c_delta_bmi.txt --qcovar $DIR/whi_white_delta_PCs.txt --out $DIR/whi_white_delta_logf2c_bmi
#gcta64 --reml --grm-bin $DIR/whi_white_delta --pheno $DIR/whi_white_delta_delta_logf2c_delta_sbp.txt --qcovar $DIR/whi_white_delta_PCs.txt --out $DIR/whi_white_delta_logf2c_sbp
#gcta64 --reml --grm-bin $DIR/whi_white_delta --pheno $DIR/whi_white_delta_delta_sfa_delta_bmi.txt --qcovar $DIR/whi_white_delta_PCs.txt --out $DIR/whi_white_delta_sfa_bmi
#gcta64 --reml --grm-bin $DIR/whi_white_delta --pheno $DIR/whi_white_delta_delta_sfa_delta_sbp.txt --qcovar $DIR/whi_white_delta_PCs.txt --out $DIR/whi_white_delta_sfa_sbp
gcta64 --reml --grm-bin $DIR/whi_white_delta --pheno $DIR/whi_white_delta_delta_sodium_delta_sbp.txt --qcovar $DIR/whi_white_delta_PCs.txt --out $DIR/whi_white_delta_sodium_sbp


### BLACKS ONLY ###
#plink2 --pfile $GENODIR/whi --keep $DIR/whi_black_delta_phenos.txt --maf 0.01 --hard-call-threshold 0 --geno 0 --make-grm-bin --out $DIR/whi_black_delta

python << EOF
import pandas as pd

phenos = pd.read_csv("$DIR/whi_black_delta_phenos.txt", sep=" ")
phenos.filter(["FID", "IID", "delta_logf2c_delta_bmi_INT"]).to_csv("$DIR/whi_black_delta_delta_logf2c_delta_bmi.txt", sep=" ", index=False)
phenos.filter(["FID", "IID", "delta_logf2c_delta_sbp_INT"]).to_csv("$DIR/whi_black_delta_delta_logf2c_delta_sbp.txt", sep=" ", index=False)
phenos.filter(["FID", "IID", "delta_sfa_delta_bmi_INT"]).to_csv("$DIR/whi_black_delta_delta_sfa_delta_bmi.txt", sep=" ", index=False)
phenos.filter(["FID", "IID", "delta_sfa_delta_sbp_INT"]).to_csv("$DIR/whi_black_delta_delta_sfa_delta_sbp.txt", sep=" ", index=False)
phenos.filter(["FID", "IID", "delta_sodium_delta_sbp_INT"]).to_csv("$DIR/whi_black_delta_delta_sodium_delta_sbp.txt", sep=" ", index=False)
phenos.filter(["FID", "IID"] + ["PC" + str(idx) for idx in range(1,6)]).to_csv("$DIR/whi_black_delta_PCs.txt", sep=" ", index=False)
EOF

#gcta64 --reml --grm-bin $DIR/whi_black_delta --pheno $DIR/whi_black_delta_delta_logf2c_delta_bmi.txt --qcovar $DIR/whi_black_delta_PCs.txt --out $DIR/whi_black_delta_logf2c_bmi
#gcta64 --reml --grm-bin $DIR/whi_black_delta --pheno $DIR/whi_black_delta_delta_logf2c_delta_sbp.txt --qcovar $DIR/whi_black_delta_PCs.txt --out $DIR/whi_black_delta_logf2c_sbp
#gcta64 --reml --grm-bin $DIR/whi_black_delta --pheno $DIR/whi_black_delta_delta_sfa_delta_bmi.txt --qcovar $DIR/whi_black_delta_PCs.txt --out $DIR/whi_black_delta_sfa_bmi
#gcta64 --reml --grm-bin $DIR/whi_black_delta --pheno $DIR/whi_black_delta_delta_sfa_delta_sbp.txt --qcovar $DIR/whi_black_delta_PCs.txt --out $DIR/whi_black_delta_sfa_sbp
gcta64 --reml --grm-bin $DIR/whi_black_delta --pheno $DIR/whi_black_delta_delta_sodium_delta_sbp.txt --qcovar $DIR/whi_black_delta_PCs.txt --out $DIR/whi_black_delta_sodium_sbp

#### ALL ETHNICITIES ###
#plink2 --pfile $GENODIR/whi --keep $DIR/whi_all_delta_phenos.txt --maf 0.01 --hard-call-threshold 0 --geno 0 --make-grm-bin --out $DIR/whi_all_delta
#
#python << EOF
#import pandas as pd
#
#phenos = pd.read_csv("$DIR/whi_all_delta_phenos.txt", sep=" ")
#phenos.filter(["FID", "IID", "delta_logf2c_delta_bmi_INT"]).to_csv("$DIR/whi_all_delta_delta_logf2c_delta_bmi.txt", sep=" ", index=False)
#phenos.filter(["FID", "IID", "delta_logf2c_delta_sbp_INT"]).to_csv("$DIR/whi_all_delta_delta_logf2c_delta_sbp.txt", sep=" ", index=False)
#phenos.filter(["FID", "IID", "delta_sfa_delta_bmi_INT"]).to_csv("$DIR/whi_all_delta_delta_sfa_delta_bmi.txt", sep=" ", index=False)
#phenos.filter(["FID", "IID", "delta_sfa_delta_sbp_INT"]).to_csv("$DIR/whi_all_delta_delta_sfa_delta_sbp.txt", sep=" ", index=False)
#phenos.filter(["FID", "IID"] + ["PC" + str(idx) for idx in range(1,6)]).to_csv("$DIR/whi_all_delta_PCs.txt", sep=" ", index=False)
#EOF
#
#gcta64 --reml --grm-bin $DIR/whi_all_delta --pheno $DIR/whi_all_delta_delta_logf2c_delta_bmi.txt --qcovar $DIR/whi_all_delta_PCs.txt --out $DIR/whi_all_delta_logf2c_bmi
#gcta64 --reml --grm-bin $DIR/whi_all_delta --pheno $DIR/whi_all_delta_delta_logf2c_delta_sbp.txt --qcovar $DIR/whi_all_delta_PCs.txt --out $DIR/whi_all_delta_logf2c_sbp
#gcta64 --reml --grm-bin $DIR/whi_all_delta --pheno $DIR/whi_all_delta_delta_sfa_delta_bmi.txt --qcovar $DIR/whi_all_delta_PCs.txt --out $DIR/whi_all_delta_sfa_bmi
#gcta64 --reml --grm-bin $DIR/whi_all_delta --pheno $DIR/whi_all_delta_delta_sfa_delta_sbp.txt --qcovar $DIR/whi_all_delta_PCs.txt --out $DIR/whi_all_delta_sfa_sbp
#
