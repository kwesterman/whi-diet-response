import sys
import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.formula.api as smf
from sklearn.preprocessing import scale
import statsmodels.formula.api as smf


# Phenotypes

phenos = (pd.read_csv(
    "../data/raw/mesa/phen/phs000209.v10.pht001116.v7.p2.c1.MESA_Exam1Main.GRU.txt",
    sep="\t", na_values=' ', skiprows=10)
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
             'pclpfn1c': 'pufa_pct', 'tcarbn1c': 'carb', 'enrgyn1c': 'tot_cal',
             'sf160n1c': 'palmitic_pct',
             'pf182n1c': 'n6',
             'pf183n1c': 'ala', 'pf205n1c': 'epa', 'pf225n1c': 'dpa',
             'pf226n1c': 'dha'}, axis=1)
    .filter(['subjID', 'sfa', 'mufa', 'pufa', 'sfa_pct', 'mufa_pct',
             'pufa_pct', 'palmitic_pct', 'tot_cal', 'carb',
             'n6', 'ala', 'epa', 'dpa', 'dha'])
    .apply(pd.to_numeric, errors="coerce")
    .assign(tot_fat = lambda x: x.sfa + x.mufa + x.pufa,
            tot_fat_pct = lambda x: x.tot_fat / x.tot_cal,
            f2c = lambda x: x.tot_fat / x.carb)
    #.assign(sfa_pct = lambda x: x.ldl.to_numeric(errors="coerce"))
    #.dropna(subset='sfa_pct')
    .assign(sfa_pct = lambda x: pd.to_numeric(x.sfa_pct, errors="coerce"),
            sfa2pufa = lambda x: x.sfa_pct / x.pufa_pct,
            n3 = lambda x: x.ala + x.epa + x.dpa + x.dha,
            n62n3 = lambda x: x.n6 / x.n3)
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

gwas_variables = ['subjID', 'ldl', 'sfa_pct', 'pufa_pct', 'f2c', 'tot_fat_pct',
                  'carb', 'palmitic_pct', 'age', 'bmi', 'tg', 'glu', 'sbp', 
                  'ldl_resid', 'sfa_pct_resid', 'race', 'sfa2pufa', 'n62n3']
phenotypes = ['sfa_ldl', 'fat_ldl', 'f2c_ldl', 'palm_ldl',
              'sfa_ldlR', 'sfaR_ldl', 'sfaR_ldlR',
              'sfa_bmi', 'f2c_bmi', 'sfa_sbp', 'f2c_sbp', 'f2c_tg', 'f2c_glu',
              'sfa2pufa_ldlR', 'sfa2pufa_tg', 'n62n3_ldlR', 'n62n3_tg']

mesa_metadata_nomeds = mesa_metadata.query('statin == False').copy()
mesa_metadata_nomeds["ldl_resid"] = smf.ols('ldl ~ age + sex + bmi', 
                                             data=mesa_metadata_nomeds).fit().resid
mesa_metadata_nomeds["sfa_pct_resid"] = smf.ols('sfa_pct ~ pufa_pct', 
                                                 data=mesa_metadata_nomeds).fit().resid

gwas_phenos = (mesa_metadata_nomeds
               .query('statin == False')
               .filter(gwas_variables + ['sex'])
               .assign(ldl = lambda x: pd.to_numeric(x.ldl, errors="coerce"))
               .dropna()
               .assign(FID = lambda x: x.subjID, IID = lambda x: x.subjID)
               .assign(sfa_ldl = lambda x: scale(x.sfa_pct) * scale(x.ldl),
                       sfa_ldlR = lambda x: scale(x.sfa_pct) *
                       scale(x.ldl_resid),
                       sfaR_ldl = lambda x: scale(x.sfa_pct_resid) *
                       scale(x.ldl),
                       sfaR_ldlR = lambda x: scale(x.sfa_pct_resid) *
                       scale(x.ldl_resid),
                       fat_ldl = lambda x: scale(x.tot_fat_pct) * scale(x.ldl),
                       f2c_ldl = lambda x: scale(x.f2c) * scale(x.ldl),
                       palm_ldl = lambda x: scale(x.palmitic_pct) *
                       scale(x.ldl),
                       sfa_bmi = lambda x: scale(x.sfa_pct) * scale(x.bmi),
                       f2c_bmi = lambda x: scale(x.f2c) * scale(x.bmi),
                       f2c_tg = lambda x: scale(x.f2c) * scale(x.tg),
                       f2c_glu = lambda x: scale(x.f2c) * scale(x.glu),
                       sfa_sbp = lambda x: scale(x.sfa_pct) * scale(x.sbp),
                       f2c_sbp = lambda x: scale(x.f2c) * scale(x.sbp),
                       sfa2pufa_ldlR = lambda x: scale(x.sfa2pufa) *
                       scale(x.ldl_resid),
                       sfa2pufa_tg = lambda x: scale(x.sfa2pufa) * scale(x.tg),
                       n62n3_ldlR = lambda x: scale(x.n62n3) *
                       scale(x.ldl_resid),
                       n62n3_tg = lambda x: scale(x.n62n3) *
                       scale(x.tg))
               .filter(['FID', 'IID', 'sex'] + phenotypes + gwas_variables))
gwas_phenos_WIN = (gwas_phenos
                   .filter(phenotypes)
                   .rename({p: p + "_WIN" for p in phenotypes}, axis=1))
gwas_phenos_INT = (gwas_phenos
                   .filter(phenotypes)
                   .apply(stats.mstats.winsorize, limits=[0.01, 0.01])
                   .apply(rank_INT)
                   .rename({p: p + "_INT" for p in phenotypes}, axis=1))
gwas_phenos_BIN = (gwas_phenos_WIN
                   .apply(lambda x: np.where(
                       x < x.quantile(0.25), 0,
                       np.where(x > x.quantile(0.75), 1, np.nan)))
                   #.apply(lambda x: np.where(~x.between(*x.quantile([0.25,
                   #                                                 0.75]))))
                   .rename({(p + "_WIN"): p + "_BIN" for p in phenotypes}, 
                           axis=1))
gwas_phenos = pd.concat([gwas_phenos, gwas_phenos_WIN, gwas_phenos_INT,
                         gwas_phenos_BIN], axis=1) 
               
gwas_phenos.query('race == "white"').to_csv(
    "../data/processed/mesa/mesa_white_gwas_phenos.txt", sep=" ", na_rep="NA", 
    index=False)
gwas_phenos.query('race == "black"').to_csv(
    "../data/processed/mesa/mesa_black_gwas_phenos.txt", sep=" ", na_rep="NA",
    index=False)
gwas_phenos.query('race == "hispanic"').to_csv(
    "../data/processed/mesa/mesa_hispanic_gwas_phenos.txt", sep=" ", 
    na_rep="NA", index=False)
gwas_phenos.query('race == "asian"').to_csv(
    "../data/processed/mesa/mesa_asian_gwas_phenos.txt", sep=" ", na_rep="NA",
    index=False)
