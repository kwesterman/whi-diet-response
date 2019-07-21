import sys
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.preprocessing import scale
import statsmodels.formula.api as smf
import itertools


# Risk factors

var_names = {'shareid': 'subjID', 'SEX': 'sex', 
             'AGE5': 'age', 'BMI5': 'bmi', 'CURRSMK5': 'smk_now', 'CPD5': 'cig_per_day',  
             'SBP5': 'sbp', 'FASTING_BG5': 'glu', 'chol': 'TC5', 'CALC_LDL5': 'ldl', 'HDL5': 'hdl', 
             'TRIG5': 'tg', 'HRX5': 'ht_med', 'LIPRX5': 'lipid_med', 'DMRX5': 'dm_med', 
             'AGE8': 'age_8', 'BMI8': 'bmi_8', 'CURRSMK8': 'smk_now_8', 'CPD8': 'cig_per_day_8',  
             'SBP8': 'sbp_8', 'FASTING_BG8': 'glu_8', 'TC8': 'chol_8', 'CALC_LDL8': 'ldl_8', 'HDL8': 'hdl_8', 
             'TRIG8': 'tg_8', 'HRX8': 'ht_med_8', 'LIPRX8': 'lipid_med_8', 'DMRX8': 'dm_med_8'} 
             
phenos_fhs_c1 = pd.read_csv("../data/raw/fhs/phen/phenos_fhs_c1.txt", sep="\t",
                            skiprows=10)
phenos_fhs_c2 = pd.read_csv("../data/raw/fhs/phen/phenos_fhs_c2.txt", sep="\t",
                            skiprows=10)
phenos_fhs = (pd.concat([phenos_fhs_c1, phenos_fhs_c2])
              .rename(var_names, axis='columns')
              .assign(sex = lambda x: np.where(x.sex == 1, "M", "F"),
                      race = "white")
              .loc[:, list(var_names.values()) + ['race']])

crp_fhs_c1 = pd.read_csv("../data/raw/fhs/phen/crp_c1.txt", sep="\t", skiprows=10)
crp_fhs_c2 = pd.read_csv("../data/raw/fhs/phen/crp_c2.txt", sep="\t", skiprows=10)
crp_fhs = (pd.concat([crp_fhs_c1, crp_fhs_c2])
           .rename({'shareid': 'subjID',
                    'crp': 'hsCRP'}, axis='columns'))

questionnaire_fhs_c1 = pd.read_csv("../data/raw/fhs/phen/exam8Data_fhs_c1.txt", 
                                   sep="\t", skiprows=10)
questionnaire_fhs_c2 = pd.read_csv("../data/raw/fhs/phen/exam8Data_fhs_c2.txt", 
                                   sep="\t", skiprows=10)
questionnaire_fhs = (pd.concat([questionnaire_fhs_c1, questionnaire_fhs_c2])
                     .rename({'shareid': 'subjID', 'H065': 'cig_start_age'},
                             axis='columns')
                     .loc[:, ['subjID', 'cig_start_age']])

pheno_data_fhs = (
    phenos_fhs 
    .merge(crp_fhs, on="subjID", how="left") 
    .merge(questionnaire_fhs, on="subjID", how="left") 
    .assign(smk_py=lambda x: np.where(x.smk_now, x.cig_per_day / 20 *
                                        (x.age - x.cig_start_age), 0),
            smk_py_8=lambda x: np.where(x.smk_now_8, x.cig_per_day_8 / 20 * 
                                        (x.age_8 - x.cig_start_age), 0))
)

# Diet

ffq_fhs_ex5_c1 = pd.read_csv("../data/raw/fhs/diet/ffq_ex5_c1.txt", sep="\t",
                             skiprows=10)
ffq_fhs_ex5_c2 = pd.read_csv("../data/raw/fhs/diet/ffq_ex5_c2.txt", sep="\t",
                             skiprows=10)
rename_diet = {'shareid': 'subjID', 'NUT_CALOR': 'tot_cal', 'NUT_SATFAT': 'sfa',
               'NUT_MONFAT': 'mufa', 'NUT_POLY': 'pufa', 'NUT_CARBO': 'carb', 
               'NUT_N3': 'n3', 'NUT_N6': 'n6', 'NUT_F182': 'linoleic', 
               'NUT_F160': 'palm'}
ffq_fhs_ex5 = (pd.concat([ffq_fhs_ex5_c1, ffq_fhs_ex5_c2])
               .rename(rename_diet, axis='columns')
               .filter(list(rename_diet.values()))
               .assign(sfa_pct = lambda x: x.sfa * 9 / x.tot_cal, 
                       mufa_pct = lambda x: x.mufa * 9 / x.tot_cal, 
                       pufa_pct = lambda x: x.pufa * 9 / x.tot_cal, 
                       sfa2pufa = lambda x: x.sfa / x.pufa,
                       palm_pct = lambda x: x.palm * 9 / 
                       x.tot_cal,
                       tot_fat = lambda x: x.sfa + x.mufa + x.pufa,
                       tot_fat_pct = lambda x: x.tot_fat / x.tot_cal,
                       f2c = lambda x: np.log(x.tot_fat / x.carb),
                       n62n3 = lambda x: x.n6 / x.n3))

ffq_fhs_ex8_c1 = pd.read_csv("../data/raw/fhs/diet/ffq_ex8_c1.txt", sep="\t",
                             skiprows=10)
ffq_fhs_ex8_c2 = pd.read_csv("../data/raw/fhs/diet/ffq_ex8_c2.txt", sep="\t",
                             skiprows=10)
ffq_fhs_ex8 = (pd.concat([ffq_fhs_ex8_c1, ffq_fhs_ex8_c2])
               .rename({'shareid': 'subjID',
                        'NUT_CALOR': 'tot_cal_8', 'NUT_SATFAT': 'sfa_8',
                        'NUT_MONFAT': 'mufa_8', 'NUT_POLY': 'pufa_8',
                        'NUT_CARBO': 'carb_8', 'NUT_PFN302': 'n3_8', 
                        'NUT_PFN602': 'n6_8',
                        'NUT_F182': 'linoleic_8', 'NUT_F160': 'palm_8'},
                       axis='columns')
               # the below could be written more elegantly
               .loc[:, lambda x: (x.columns == 'subjID') | 
                    x.columns.str.contains('_8')]
               .assign(sfa_pct_8 = lambda x: x.sfa_8 * 9 / x.tot_cal_8, 
                       mufa_pct_8 = lambda x: x.mufa_8 * 9 / x.tot_cal_8, 
                       pufa_pct_8 = lambda x: x.pufa_8 * 9 / x.tot_cal_8, 
                       sfa2pufa_8 = lambda x: x.sfa_8 / x.pufa_8,
                       palm_pct_8 = lambda x: x.palm_8 * 9 / 
                       x.tot_cal_8,
                       tot_fat_8 = lambda x: x.sfa_8 + x.mufa_8 + x.pufa_8,
                       tot_fat_pct_8 = lambda x: x.tot_fat_8 / x.tot_cal_8,
                       f2c_8 = lambda x: np.log(x.tot_fat_8 / x.carb_8),
                       n62n3_8 = lambda x: x.n6_8 / x.n3_8))

ffq_fhs = pd.merge(ffq_fhs_ex5, ffq_fhs_ex8, on="subjID", how="outer")

# Merge and save
               
fhs_metadata = pd.merge(pheno_data_fhs, ffq_fhs, on="subjID")
fhs_metadata.to_csv("../data/processed/metadata_fhs.csv", index=False)

# GWAS prep

# Pedigree
ped = (pd.read_csv("../data/raw/fhs/fhs_pedigree.txt", sep="\t", skiprows=10)
       .loc[:, ['pedno', 'shareid']]
       .rename(columns={'pedno': 'FID', 'shareid': 'IID'}))

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

fhs_metadata_medsAdj = (fhs_metadata
                        .assign(ldl = lambda x: np.where(x.lipid_med, x.ldl /
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
    fhs_metadata_medsAdj[rf + "_R"] = smf.ols(
        rf + ' ~ age + sex + bmi + race',
        data=fhs_metadata_medsAdj).fit().resid
risk_factors.extend([rf + "_R" for rf in risk_factors])

gwas_covars = ['subjID', 'sex', 'age', 'bmi']
gwas_covars.extend(diet_vars + risk_factors)

dv_rf_combos = list(itertools.product(diet_vars, risk_factors))
dv_rf_combo_names = [dv + "_" + rf for dv, rf in dv_rf_combos]

gwas_phenos = (fhs_metadata_medsAdj
               #.query('(lipid_med == False | lipid_med.isna()) & dm_trial == False')
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

gwas_phenos.to_csv("../data/processed/fhs/fhs_gwas_phenos_statinAdj.txt", sep=" ", 
                   na_rep="NA", index=False)
