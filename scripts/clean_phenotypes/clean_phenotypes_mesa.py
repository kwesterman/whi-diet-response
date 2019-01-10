import sys
import numpy as np
import pandas as pd
from scipy import stats
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
             'pclpfn1c': 'pufa_pct', 'tcarbn1c': 'carb', 'enrgyn1c': 'tot_cal'}, axis=1)
    .filter(['subjID', 'sfa', 'mufa', 'pufa', 'sfa_pct', 'mufa_pct',
             'pufa_pct', 'tot_cal', 'carb', 'f2c'])
    .apply(pd.to_numeric, errors="coerce")
    .assign(tot_fat = lambda x: x.sfa + x.mufa + x.pufa,
            tot_fat_pct = lambda x: x.tot_fat / x.tot_cal,
            f2c = lambda x: x.tot_fat / x.carb)
    #.assign(sfa_pct = lambda x: x.ldl.to_numeric(errors="coerce"))
    #.dropna(subset='sfa_pct')
    .assign(sfa_pct = lambda x: pd.to_numeric(x.sfa_pct, errors="coerce"))
    .dropna(subset=["sfa_pct"]))

# Merge and save

mesa_metadata = pd.merge(phenos, ffq, on="subjID")
mesa_metadata.to_csv("../data/processed/metadata_mesa.csv", index=False)


# GWAS prep

# Construct outcome feature and prepare plink-friendly dataset
def rank_INT(values):
    # Computes rank-based inverse normal transformation of a Pandas series
    c = 3.0 / 8
    ranks = stats.rankdata(values)
    product_INT = stats.norm.ppf((ranks - c) / (len(ranks) - 2*c + 1))
    return product_INT

def winsorize(values, SDs_from_mean=4):
    # Winsorize values based on SDs from mean (rather than percentiles)
    sd = np.std(values)
    values[values < values.mean() - 4 * sd] = values.mean() - 4 * sd
    values[values > values.mean() + 4 * sd] = values.mean() + 4 * sd
    return values

gwas_variables = ['subjID', 'ldl', 'sfa_pct', 'pufa_pct', 'f2c', 'tot_fat_pct',
                  'carb', 'age', 'bmi', 'tg', 'glu', 'race']
phenotypes = ['sfa_ldl', 'fat_ldl', 'f2c_ldl', 
              'sfa_bmi', 'f2c_bmi', 'f2c_tg', 'f2c_glu']
gwas_phenos = (mesa_metadata
               .query('statin == False')
               .filter(gwas_features + ['sex'])
               .assign(ldl = lambda x: pd.to_numeric(x.ldl, errors="coerce"))
               .dropna()
               .assign(FID = lambda x: x.subjID, IID = lambda x: x.subjID)
               .assign(sfa_ldl = lambda x: scale(x.sfa_pct) * scale(x.ldl),
                       fat_ldl = lambda x: scale(x.tot_fat_pct) * scale(x.ldl),
                       f2c_ldl = lambda x: scale(x.f2c) * scale(x.ldl),
                       sfa_bmi = lambda x: scale(x.sfa_pct) * scale(x.bmi),
                       f2c_bmi = lambda x: scale(x.f2c) * scale(x.bmi),
                       f2c_tg = lambda x: scale(x.f2c) * scale(x.tg),
                       f2c_glu = lambda x: scale(x.f2c) * scale(x.glu))
               .filter(['FID', 'IID', 'sex'] + phenotypes + gwas_features))
gwas_phenos_WIN = (gwas_phenos
                   .filter(phenotypes)
                   .apply(winsorize)
                   .rename({p: p + "_WIN" for p in phenotypes}, axis=1))
gwas_phenos_INT = (gwas_phenos
                   .filter(phenotypes)
                   .apply(rank_INT)
                   .rename({p: p + "_INT" for p in phenotypes}, axis=1))
gwas_phenos = pd.concat([gwas_phenos, gwas_phenos_WIN, gwas_phenos_INT],
                        axis=1) 
               
gwas_phenos.query('race == "white"').to_csv(
    "../data/processed/mesa/mesa_white_gwas_phenos.txt", sep=" ", index=False)
gwas_phenos.query('race == "black"').to_csv(
    "../data/processed/mesa/mesa_black_gwas_phenos.txt", sep=" ", index=False)
gwas_phenos.query('race == "hispanic"').to_csv(
    "../data/processed/mesa/mesa_hispanic_gwas_phenos.txt", sep=" ", index=False)
gwas_phenos.query('race == "asian"').to_csv(
    "../data/processed/mesa/mesa_asian_gwas_phenos.txt", sep=" ", index=False)
