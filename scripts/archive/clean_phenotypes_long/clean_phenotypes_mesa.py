import sys
import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.formula.api as smf
from sklearn.preprocessing import scale
import statsmodels.formula.api as smf
import itertools


# Phenotypes

phenos = (pd.read_csv( "../data/raw/mesa/phen/SHARe_Exam1Main.txt", sep="\t", na_values=' ')
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

phenos_5 = (pd.read_csv( "../data/raw/mesa/phen/SHARe_Exam5Main.txt", sep="\t", na_values=' ')
    .rename({'sidno': 'subjID', 'age1c': 'age', 'gender1': 'sex', 
             'race5c': 'race', 'bmi5c': 'bmi', 'cig5c': 'smk_now', 
             'pkyrs5c': 'smk_py', 
             'dm031t': 'diabetes_recode', 'diabhx5': 'dm_med',
             'htn5c': 'hypertension', 'htnmed5c': 'ht_med', 'ldl5': 'ldl', 'hdl5': 'hdl',
             'chol5': 'chol', 'trig5': 'tg', 'glucose5': 'glu', 'sbp5c': 'sbp',
             'sttn5c': 'statin'}, axis=1)
    .assign(sex = lambda x: np.where(x.sex == 0, "F", "M"),
            statin = lambda x: pd.to_numeric(x.statin, errors="coerce"))
    .replace({'race': {1: 'white', 2: 'asian', 3: 'black', 4: 'hispanic'}})
    .filter(['subjID', 'sex', 'age', 'race', 'bmi', 'smk_now', 'smk_py',
             'diabetes_recode', 'dm_med', 'hypertension', 'ht_med', 'statin', 
             'ldl', 'hdl', 'chol', 'tg', 'glu', 'sbp', 'hsCRP']))

phenos_with_deltas = (pd.merge(phenos, phenos_5, on="subjID", 
                               suffixes=("", "_5"))
                      .assign(delta_chol=lambda x: x.chol_5 - x.chol,
                              delta_ldl=lambda x: x.ldl_5 - x.ldl,
                              delta_hdl=lambda x: x.hdl_5 - x.hdl,
                              delta_tg=lambda x: x.tg_5 - x.tg,
                              delta_glu=lambda x: x.glu_5 - x.glu,
                              delta_sbp=lambda x: x.sbp_5 - x.sbp,
                              delta_bmi=lambda x: x.bmi_5 - x.bmi))


# FFQ

def winsorize(x, thresh):
    return x.clip(lower=x.quantile(thresh), upper=x.quantile(1 - thresh))

ffq_nut = (pd.read_csv("../data/raw/mesa/diet/SHARe_Exam1DietNutrients.txt",
                   sep="\t")
       .rename({'sidno': 'subjID', 'tsfan1c': 'sfa', 'tmufan1c': 'mufa',
                'tfatn1c': 'fat', 'tprtnn1c': 'pro', 'tfibrn1c': 'fiber',
                'sf160n1c': 'palm', 'alcn1c': 'alc',
             'tpufan1c': 'pufa', 'pclsfn1c': 'sfa_pct', 'pclmfn1c': 'mufa_pct',
             'pclpfn1c': 'pufa_pct', 'tcarbn1c': 'cho', 'enrgyn1c': 'tot_cal',
             #'sf160n1c': 'palm_pct',
             'pf182n1c': 'n6',
             'pf183n1c': 'ala', 'pf205n1c': 'epa', 'pf225n1c': 'dpa',
             'pf226n1c': 'dha'}, axis=1)
    .filter(['subjID', 'sfa', 'mufa', 'pufa', 'palm', 'fat', 'cho', 'fiber',
             'sfa_pct', 'mufa_pct', 'pro', 'alc',
             'pufa_pct', 'palm_pct', 'tot_cal', 
             'n6', 'ala', 'epa', 'dpa', 'dha'])
    .apply(pd.to_numeric, errors="coerce")
    #.assign(sfa = lambda x: x.sfa_pct,
    #        mufa = lambda x: x.mufa_pct,
    #        pufa = lambda x: x.pufa_pct,
    #        tot_fat = lambda x: x.sfa + x.mufa + x.pufa,
    #        palm = lambda x: x.palm_pct,
    .assign(logf2c = lambda x: np.log(x.fat / x.cho))
    .assign(sfa = lambda x: pd.to_numeric(x.sfa, errors="coerce"),
            logsfa2pufa = lambda x: np.log(x.sfa / x.pufa),
            n3 = lambda x: x.ala + x.epa + x.dpa + x.dha,
            logn62n3 = lambda x: np.log(x.n6 / x.n3))
    .dropna(subset=["sfa"]))

ffq_items = (pd.read_csv("../data/raw/mesa/diet/phs000209.v10.pht001115.v3.p2.c1.MESA_Exam1Diet.GRU.txt", 
                         sep="\t", skiprows=10)
             .rename({'sidno': 'subjID', 'fgvgreenleafy1c': 'GLVEG',
                      'fghfprocmeat1c': 'PROCMEAT', 'fgwholegrain1c': 'WG',
                      'fgdesserts1c': 'DESSERT', 'fgseedsnuts1c':
                      'seeds_nuts'}, axis=1)
             .filter(['subjID', 'GLVEG', 'PROCMEAT', 'WG', 'DESSERT',
                      'seeds_nuts'])
             .apply(pd.to_numeric, errors="coerce")
             .dropna()
             .assign(GLVEG=lambda x: scale(winsorize(x.GLVEG, 0.05)),
                     PROCMEAT=lambda x: scale(winsorize(x.PROCMEAT, 0.05)),
                     WG=lambda x: scale(winsorize(x.WG, 0.05)),
                     DESSERT=lambda x: scale(winsorize(x.DESSERT, 0.05)),
                     seeds_nuts=lambda x: scale(winsorize(x.seeds_nuts, 0.05)),
                     myhei=lambda x: 2 * x.GLVEG + x.WG + x.seeds_nuts -
                     x.PROCMEAT - x.DESSERT))

ffq = pd.merge(ffq_nut, ffq_items, on="subjID")

ffq_nut_5 = (pd.read_csv("../data/raw/mesa/diet/SHARe_Exam5DietNutrients.txt",
                   sep="\t")
       .rename({'sidno': 'subjID', 'tsfan5c': 'sfa', 'tmufan5c': 'mufa',
             'tpufan5c': 'pufa', 'pclsfn5c': 'sfa_pct', 'pclmfn5c': 'mufa_pct',
             'pclpfn5c': 'pufa_pct', 'tcarbn5c': 'cho', 'enrgyn5c': 'tot_cal',
             'tfatn5c': 'fat', 'tprtnn5c': 'pro', 'tfibrn5c': 'fiber',
                'sf160n5c': 'palm', 'alcn5c': 'alc',
             #'sf160n5c': 'palm_pct',
             'pf182n5c': 'n6',
             'pf183n5c': 'ala', 'pf205n5c': 'epa', 'pf225n5c': 'dpa',
             'pf226n5c': 'dha'}, axis=1)
    .filter(['subjID', 'sfa', 'mufa', 'pufa', 'palm', 'fat', 'cho', 'fiber',
             'sfa_pct', 'mufa_pct', 'pro', 'alc',
             'pufa_pct', 'palm_pct', 'tot_cal', 
             'n6', 'ala', 'epa', 'dpa', 'dha'])
    .apply(pd.to_numeric, errors="coerce")
             .assign(logf2c = lambda x: np.log(x.fat / x.cho),
                     logsfa2pufa = lambda x: np.log(x.sfa / x.pufa),
            n3 = lambda x: x.ala + x.epa + x.dpa + x.dha,
            fat = lambda x: x.sfa + x.mufa + x.pufa,
            logn62n3 = lambda x: np.log(x.n6 / x.n3))
    .dropna(subset=["sfa"]))

ffq_5 = ffq_nut_5

ffq_with_deltas = (pd.merge(ffq, ffq_5, on="subjID", how="left", suffixes=("", "_5"))
                   .assign(delta_sfa = lambda x: x.sfa_5 - x.sfa,
                           delta_pufa = lambda x: x.pufa_5 - x.pufa,
                           delta_logf2c = lambda x: x.logf2c_5 - x.logf2c,
                           delta_palm = lambda x: x.palm_5 - x.palm,
                           delta_logsfa2pufa = lambda x: x.logsfa2pufa_5 - x.logsfa2pufa,
                           delta_fat = lambda x: x.fat_5 - x.fat,
                           delta_logn62n3 = lambda x: x.logn62n3_5 - x.logn62n3))

# Merge and save

mesa_metadata = pd.merge(phenos_with_deltas, ffq_with_deltas, on="subjID")
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

diet_vars = ['sfa', 'pufa', 'fat', 'logf2c', 'palm',
             'logsfa2pufa', 'logn62n3', 'n3', 'myhei']
diet_vars.extend(["delta_" + dv for dv in diet_vars if dv not in ["n3",
                                                                  "myhei"]])
risk_factors = ["chol", "ldl", "hdl", "tg", "glu", "hsCRP", "sbp", "bmi"]
risk_factors.extend(["delta_" + rf for rf in risk_factors 
                     if rf != "hsCRP"])

#for rf in risk_factors:
#    mesa_metadata_medsAdj[rf + "_R"] = smf.ols(
#        rf + ' ~ age + sex + bmi + race',
#        data=mesa_metadata_medsAdj).fit().resid
#risk_factors.extend([rf + "_R" for rf in risk_factors])

gwas_covars = ['subjID', 'sex', 'race', 'age']
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

final_phenos = dv_rf_combo_names + diet_vars + risk_factors
gwas_phenos_INT = (gwas_phenos  # Inverse-normal transforms of phenotypes
                   .filter(final_phenos)
                   .apply(rank_INT)
                   .rename({p: p + "_INT" for p in final_phenos}, axis=1))
gwas_phenos = pd.concat([gwas_phenos, gwas_phenos_INT], axis=1) 
               
gwas_phenos.query('race == "white"').to_csv(
    "../data/processed/mesa/mesa_white_gwas_phenos_long.txt", sep=" ", na_rep="NA", 
    index=False)
gwas_phenos.query('race == "black"').to_csv(
    "../data/processed/mesa/mesa_black_gwas_phenos_long.txt", sep=" ", na_rep="NA",
    index=False)
gwas_phenos.query('race == "hispanic"').to_csv(
    "../data/processed/mesa/mesa_hispanic_gwas_phenos_long.txt", sep=" ", 
    na_rep="NA", index=False)
gwas_phenos.query('race == "asian"').to_csv(
    "../data/processed/mesa/mesa_asian_gwas_phenos_long.txt", sep=" ", na_rep="NA",
    index=False)
