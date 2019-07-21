import sys
import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.formula.api as smf
from sklearn.preprocessing import scale
import statsmodels.formula.api as smf
import itertools


# Phenotypes

phenos = (pd.read_csv(
    "../data/raw/mesa/phen/phs000209.v10.pht001116.v7.p2.c1.MESA_Exam1Main.GRU.txt",
    sep="\t", na_values=' ', skiprows=10)
    .rename({'sidno': 'subjID', 'age1c': 'age', 'gender1': 'sex', 
             'race1c': 'race', 'bmi1c': 'bmi', 'cig1c': 'smk_now', 
             'pkyrs1c': 'smk_py', 
             'dm031t': 'diabetes_recode', 'diabhx1': 'dm_med',
             'htn1c': 'hypertension', 'htnmed1c': 'ht_med', 'ldl1': 'ldl', 'hdl1': 'hdl',
             'chol1': 'chol', 'trig1': 'tg', 'glucos1c': 'glu', 'sbp1c': 'sbp',
             'crp1': 'hsCRP', 'sttn1c': 'statin'}, axis=1)
    .assign(sex = lambda x: np.where(x.sex == 0, "F", "M"),
            statin = lambda x: pd.to_numeric(x.statin, errors="coerce"))
    .replace({'race': {1: 'white', 2: 'asian', 3: 'black', 4: 'hispanic'}})
    .filter(['subjID', 'sex', 'age', 'race', 'bmi', 'smk_now', 'smk_py',
             'diabetes_recode', 'dm_med', 'hypertension', 'ht_med', 'statin', 
             'ldl', 'hdl', 'chol', 'tg', 'glu', 'sbp', 'hsCRP']))

# FFQ

ffq = (pd.read_csv(
    "../data/raw/mesa/diet/phs000209.v10.pht002107.v3.p2.c1.MESA_Exam1DietNutrients.GRU.txt",
    sep="\t", skiprows=10)
    .rename({'sidno': 'subjID', 'tsfan1c': 'sfa', 'tmufan1c': 'mufa',
             'tpufan1c': 'pufa', 'pclsfn1c': 'sfa_pct', 'pclmfn1c': 'mufa_pct',
             'pclpfn1c': 'pufa_pct', 'tcarbn1c': 'carb', 'enrgyn1c': 'tot_cal',
             'sf160n1c': 'palm_pct',
             'pf182n1c': 'n6',
             'pf183n1c': 'ala', 'pf205n1c': 'epa', 'pf225n1c': 'dpa',
             'pf226n1c': 'dha'}, axis=1)
    .filter(['subjID', 'sfa', 'mufa', 'pufa', 'sfa_pct', 'mufa_pct',
             'pufa_pct', 'palm_pct', 'tot_cal', 'carb',
             'n6', 'ala', 'epa', 'dpa', 'dha'])
    .apply(pd.to_numeric, errors="coerce")
    .assign(sfa = lambda x: x.sfa_pct,
            pufa = lambda x: x.pufa_pct,
            tot_fat = lambda x: x.sfa_pct + x.mufa_pct + x.pufa_pct,
            palm = lambda x: x.palm_pct,
            f2c = lambda x: np.log(x.tot_fat / x.carb))
    #.assign(sfa_pct = lambda x: x.ldl.to_numeric(errors="coerce"))
    #.dropna(subset='sfa_pct')
    .assign(sfa = lambda x: pd.to_numeric(x.sfa, errors="coerce"),
            sfa2pufa = lambda x: x.sfa / x.pufa,
            n3 = lambda x: x.ala + x.epa + x.dpa + x.dha,
            n62n3 = lambda x: x.n6 / x.n3)
    .dropna(subset=["sfa"]))

# Merge and save

mesa_metadata = pd.merge(phenos, ffq, on="subjID")
mesa_metadata.to_csv("../data/processed/metadata_mesa.csv", index=False)


# GWAS prep

# Construct outcome feature and prepare plink-friendly dataset
def rank_INT(values):
    # Computes rank-based inverse normal transformation of a Pandas series
    # while preserving missing values
    c = 3.0 / 8
    ranks = values.dropna().rank()
    product_INT = stats.norm.ppf((ranks - c) / (len(ranks) - 2*c + 1))
    values[~values.isna()] = product_INT
    return values

def scale(arr):
    return (arr - np.nanmean(arr)) / np.nanstd(arr)

mesa_metadata_medsAdj = (mesa_metadata
                        .assign(ldl = lambda x: np.where(x.statin, x.ldl /
                                                         0.75, x.ldl),
                                glu = lambda x: np.where(x.dm_med, x.glu /
                                                         0.75, x.glu),
                                sbp = lambda x: np.where(x.sbp, x.sbp  + 15,
                                                         x.sbp))
                        .copy())

diet_vars = ['sfa', 'pufa', 'tot_fat', 'f2c', 'palm',
             'sfa2pufa', 'n62n3']
risk_factors = ["ldl", "tg", "glu", "hsCRP", "sbp"]

for rf in risk_factors:
    mesa_metadata_medsAdj[rf + "_R"] = smf.ols(
        rf + ' ~ age + sex + bmi + race',
        data=mesa_metadata_medsAdj).fit().resid
risk_factors.extend([rf + "_R" for rf in risk_factors])

gwas_covars = ['subjID', 'sex', 'race', 'age', 'bmi']
gwas_covars.extend(diet_vars + risk_factors)

dv_rf_combos = list(itertools.product(diet_vars, risk_factors))
dv_rf_combo_names = [dv + "_" + rf for dv, rf in dv_rf_combos]

gwas_phenos = (mesa_metadata_medsAdj
               .filter(gwas_covars)
               #.dropna()
               .assign(FID = lambda x: x.subjID, IID = lambda x: x.subjID)
               .filter(['FID', 'IID', 'sex'] + dv_rf_combo_names + 
                       [v for v in gwas_covars if v != "sex"]))
for dv, rf in dv_rf_combos:  # Calculate scaled products for each diet/RF combo
    gwas_phenos[dv + "_" + rf] = scale(gwas_phenos[dv]) * scale(gwas_phenos[rf])

final_phenos = dv_rf_combo_names + ["f2c"]
gwas_phenos_INT = (gwas_phenos  # Inverse-normal transforms of phenotypes
                   .filter(final_phenos)
                   .apply(rank_INT)
                   .rename({p: p + "_INT" for p in final_phenos}, axis=1))
gwas_phenos = pd.concat([gwas_phenos, gwas_phenos_INT], axis=1) 
               
gwas_phenos.query('race == "white"').to_csv(
    "../data/processed/mesa/mesa_white_gwas_phenos_statinAdj.txt", sep=" ", na_rep="NA", 
    index=False)
gwas_phenos.query('race == "black"').to_csv(
    "../data/processed/mesa/mesa_black_gwas_phenos_statinAdj.txt", sep=" ", na_rep="NA",
    index=False)
gwas_phenos.query('race == "hispanic"').to_csv(
    "../data/processed/mesa/mesa_hispanic_gwas_phenos_statinAdj.txt", sep=" ", 
    na_rep="NA", index=False)
gwas_phenos.query('race == "asian"').to_csv(
    "../data/processed/mesa/mesa_asian_gwas_phenos_statinAdj.txt", sep=" ", na_rep="NA",
    index=False)
