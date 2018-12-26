import sys
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from sklearn.preprocessing import scale


# Phenotypes

phenos = (pd.read_csv(
    "../data/raw/mesa/phen/phs000209.v10.pht001116.v7.p2.c1.MESA_Exam1Main.GRU.txt",
    sep="\t", skiprows=10)
    .rename({'sidno': 'subjID', 'age1c': 'age', 'gender1': 'sex', 
             'race1c': 'race', 'bmi1c': 'bmi', 'cig1c': 'smk_now', 
             'pkyrs1c': 'smk_py', 'dm031t': 'diabetes_recode',
             'htn1c': 'hypertension', 'ldl1': 'ldl', 'hdl1': 'hdl',
             'chol1': 'chol', 'trig1': 'tg', 'glucos1c': 'glu', 'sbp1c': 'sbp',
             'crp1': 'crp', 'sttn1c': 'statin'}, axis=1)
    .assign(sex = lambda x: np.where(x.sex == 0, "F", "M"),
            statin = lambda x: pd.to_numeric(x.statin, errors="coerce"))
    .replace({'race': {1: 'white', 2: 'asian', 3: 'black', 4: 'hispanic'}})
    .filter(['subjID', 'sex', 'age', 'race', 'bmi', 'smk_now', 'smk_py',
             'diabetes_recode', 'hypertension', 'statin', 'ldl', 'hdl', 'chol',
             'tg', 'glu', 'sbp', 'crp']))

# FFQ

ffq = (pd.read_csv(
    "../data/raw/mesa/diet/phs000209.v10.pht002107.v3.p2.c1.MESA_Exam1DietNutrients.GRU.txt",
    sep="\t", skiprows=10)
    .rename({'sidno': 'subjID', 'tsfan1c': 'sfa', 'tmufan1c': 'mufa',
             'tpufan1c': 'pufa', 'pclsfn1c': 'sfa_pct', 'pclmfn1c': 'mufa_pct',
             'pclpfn1c': 'pufa_pct'}, axis=1)
    #.assign(sfa_pct = lambda x: x.ldl.to_numeric(errors="coerce"))
    #.dropna(subset='sfa_pct')
    .filter(['subjID', 'sfa', 'mufa', 'pufa', 'sfa_pct', 'mufa_pct',
             'pufa_pct'])
    .assign(sfa_pct = lambda x: pd.to_numeric(x.sfa_pct, errors="coerce"))
    .dropna(subset=["sfa_pct"]))

# Merge and save

mesa_metadata = pd.merge(phenos, ffq, on="subjID")
mesa_metadata.to_csv("../data/processed/metadata_mesa.csv", index=False)


# GWAS prep

# Construct outcome feature and prepare plink-friendly dataset
gwas_features = ['subjID', 'ldl', 'sfa_pct', 'pufa_pct', 'age', 'bmi', 'race']
gwas_phenos = (mesa_metadata
               .query('statin == False')
               .filter(gwas_features + ['sex'])
               .assign(ldl = lambda x: pd.to_numeric(x.ldl, errors="coerce"))
               .dropna()
               .assign(sfa_ldl_product = lambda x: scale(x.sfa_pct) *
                       scale(x.ldl))
               .assign(FID = lambda x: x.subjID,
                       IID = lambda x: x.subjID)
               .filter(['FID', 'IID', 'sex', 'sfa_ldl_product'] + gwas_features))

gwas_phenos.query('race == "white"').to_csv(
    "../data/processed/mesa/mesa_white_gwas_phenos.txt", sep=" ", index=False)
gwas_phenos.query('race == "black"').to_csv(
    "../data/processed/mesa/mesa_black_gwas_phenos.txt", sep=" ", index=False)
gwas_phenos.query('race == "hispanic"').to_csv(
    "../data/processed/mesa/mesa_hispanic_gwas_phenos.txt", sep=" ", index=False)
gwas_phenos.query('race == "asian"').to_csv(
    "../data/processed/mesa/mesa_asian_gwas_phenos.txt", sep=" ", index=False)
