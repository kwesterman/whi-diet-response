#!/bin/bash


module load plink2
module load gcta


GENODIR=../data/processed/whi
DIR=../data/processed/architecture

### WHITES ONLY ###
#plink2 --pfile $GENODIR/whi --keep $DIR/whi_white_cs_biochemical_phenos.txt --maf 0.01 --hard-call-threshold 0 --geno 0 --make-grm-bin --out $DIR/whi_white_cs_biochemical

python << EOF
import pandas as pd

phenos = pd.read_csv("$DIR/whi_white_cs_biochemical_phenos.txt", sep=" ")
phenos.filter(["FID", "IID", "logf2c_ldl_INT"]).to_csv("$DIR/whi_white_cs_logf2c_ldl.txt", sep=" ", index=False)
phenos.filter(["FID", "IID", "logf2c_hdl_INT"]).to_csv("$DIR/whi_white_cs_logf2c_hdl.txt", sep=" ", index=False)
phenos.filter(["FID", "IID", "logf2c_tg_INT"]).to_csv("$DIR/whi_white_cs_logf2c_tg.txt", sep=" ", index=False)
phenos.filter(["FID", "IID", "logf2c_glu_INT"]).to_csv("$DIR/whi_white_cs_logf2c_glu.txt", sep=" ", index=False)
phenos.filter(["FID", "IID", "logf2c_hsCRP_INT"]).to_csv("$DIR/whi_white_cs_logf2c_hsCRP.txt", sep=" ", index=False)
phenos.filter(["FID", "IID", "sfa_ldl_INT"]).to_csv("$DIR/whi_white_cs_sfa_ldl.txt", sep=" ", index=False)
phenos.filter(["FID", "IID", "sfa_hdl_INT"]).to_csv("$DIR/whi_white_cs_sfa_hdl.txt", sep=" ", index=False)
phenos.filter(["FID", "IID", "sfa_tg_INT"]).to_csv("$DIR/whi_white_cs_sfa_tg.txt", sep=" ", index=False)
phenos.filter(["FID", "IID", "sfa_glu_INT"]).to_csv("$DIR/whi_white_cs_sfa_glu.txt", sep=" ", index=False)
phenos.filter(["FID", "IID", "sfa_hsCRP_INT"]).to_csv("$DIR/whi_white_cs_sfa_hsCRP.txt", sep=" ", index=False)
phenos.filter(["FID", "IID", "n3_tg_INT"]).to_csv("$DIR/whi_white_cs_n3_tg.txt", sep=" ", index=False)
phenos.filter(["FID", "IID"] + ["PC" + str(idx) for idx in range(1,6)]).to_csv("$DIR/whi_white_cs_biochemical_PCs.txt", sep=" ", index=False)
EOF

#gcta64 --reml --grm-bin $DIR/whi_white_cs_biochemical --pheno $DIR/whi_white_cs_logf2c_ldl.txt --qcovar $DIR/whi_white_cs_biochemical_PCs.txt --out $DIR/whi_white_cs_logf2c_ldl
#gcta64 --reml --grm-bin $DIR/whi_white_cs_biochemical --pheno $DIR/whi_white_cs_logf2c_hdl.txt --qcovar $DIR/whi_white_cs_biochemical_PCs.txt --out $DIR/whi_white_cs_logf2c_hdl
#gcta64 --reml --grm-bin $DIR/whi_white_cs_biochemical --pheno $DIR/whi_white_cs_logf2c_tg.txt --qcovar $DIR/whi_white_cs_biochemical_PCs.txt --out $DIR/whi_white_cs_logf2c_tg
#gcta64 --reml --grm-bin $DIR/whi_white_cs_biochemical --pheno $DIR/whi_white_cs_logf2c_glu.txt --qcovar $DIR/whi_white_cs_biochemical_PCs.txt --out $DIR/whi_white_cs_logf2c_glu
#gcta64 --reml --grm-bin $DIR/whi_white_cs_biochemical --pheno $DIR/whi_white_cs_logf2c_hsCRP.txt --qcovar $DIR/whi_white_cs_biochemical_PCs.txt --out $DIR/whi_white_cs_logf2c_hsCRP
#gcta64 --reml --grm-bin $DIR/whi_white_cs_biochemical --pheno $DIR/whi_white_cs_sfa_ldl.txt --qcovar $DIR/whi_white_cs_biochemical_PCs.txt --out $DIR/whi_white_cs_sfa_ldl
#gcta64 --reml --grm-bin $DIR/whi_white_cs_biochemical --pheno $DIR/whi_white_cs_sfa_hdl.txt --qcovar $DIR/whi_white_cs_biochemical_PCs.txt --out $DIR/whi_white_cs_sfa_hdl
#gcta64 --reml --grm-bin $DIR/whi_white_cs_biochemical --pheno $DIR/whi_white_cs_sfa_tg.txt --qcovar $DIR/whi_white_cs_biochemical_PCs.txt --out $DIR/whi_white_cs_sfa_tg
#gcta64 --reml --grm-bin $DIR/whi_white_cs_biochemical --pheno $DIR/whi_white_cs_sfa_glu.txt --qcovar $DIR/whi_white_cs_biochemical_PCs.txt --out $DIR/whi_white_cs_sfa_glu
#gcta64 --reml --grm-bin $DIR/whi_white_cs_biochemical --pheno $DIR/whi_white_cs_sfa_hsCRP.txt --qcovar $DIR/whi_white_cs_biochemical_PCs.txt --out $DIR/whi_white_cs_sfa_hsCRP
gcta64 --reml --grm-bin $DIR/whi_white_cs_biochemical --pheno $DIR/whi_white_cs_n3_tg.txt --qcovar $DIR/whi_white_cs_biochemical_PCs.txt --out $DIR/whi_white_cs_n3_tg
