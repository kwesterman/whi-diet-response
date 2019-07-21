import sys
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.preprocessing import scale
import statsmodels.formula.api as smf


# Risk factors

var_names = {'shareid': 'subjID', 'SEX': 'sex', 
             'AGE5': 'age_5', 'BMI5': 'bmi_5', 'CURRSMK5': 'smk_now_5', 'CPD5': 'cig_per_day_5',  
             'SBP5': 'sbp_5', 'FASTING_BG5': 'glu_5', 'chol_5': 'TC5', 'CALC_LDL5': 'ldl_5', 'HDL5': 'hdl_5', 
             'TRIG5': 'tg_5', 'HRX5': 'ht_med_5', 'LIPRX5': 'lipid_med_5', 'DMRX5': 'dm_med_5', 
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
           .rename({'shareid': 'subjID'}, axis='columns'))

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
    .merge(crp_fhs, on="subjID") 
    .merge(questionnaire_fhs, on="subjID") 
    .assign(smk_py_5=lambda x: np.where(x.smk_now_5, x.cig_per_day_5 / 20 *
                                        (x.age_5 - x.cig_start_age), 0),
            smk_py_8=lambda x: np.where(x.smk_now_8, x.cig_per_day_8 / 20 * 
                                        (x.age_8 - x.cig_start_age), 0))
)

# Diet

ffq_fhs_ex5_c1 = pd.read_csv("../data/raw/fhs/diet/ffq_ex5_c1.txt", sep="\t",
                             skiprows=10)
ffq_fhs_ex5_c2 = pd.read_csv("../data/raw/fhs/diet/ffq_ex5_c2.txt", sep="\t",
                             skiprows=10)
ffq_fhs_ex5 = (pd.concat([ffq_fhs_ex5_c1, ffq_fhs_ex5_c2])
               .rename({'shareid': 'subjID',
                        'NUT_CALOR': 'tot_cal_5', 'NUT_SATFAT': 'sfa_5',
                        'NUT_MONFAT': 'mufa_5', 'NUT_POLY': 'pufa_5',
                        'NUT_CARBO': 'carb_5', 'NUT_N3': 'n3_5', 'NUT_N6':
                        'n6_5',
                        'NUT_F182': 'linoleic_5', 'NUT_F160': 'palmitic_5'},
                       axis='columns')
               .iloc[:, lambda x: (x.columns == 'subjID') | 
                     x.columns.str.contains('_5')]
               .assign(sfa_pct_5 = lambda x: x.sfa_5 * 9 / x.tot_cal_5, 
                       mufa_pct_5 = lambda x: x.mufa_5 * 9 / x.tot_cal_5, 
                       pufa_pct_5 = lambda x: x.pufa_5 * 9 / x.tot_cal_5, 
                       sfa2pufa_5 = lambda x: x.sfa_5 / x.pufa_5,
                       palmitic_pct_5 = lambda x: x.palmitic_5 * 9 / 
                       x.tot_cal_5,
                       tot_fat_5 = lambda x: x.sfa_5 + x.mufa_5 + x.pufa_5,
                       tot_fat_pct_5 = lambda x: x.tot_fat_5 / x.tot_cal_5,
                       f2c_5 = lambda x: x.tot_fat_5 / x.carb_5,
                       n62n3_5 = lambda x: x.n6_5 / x.n3_5))

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
                        'NUT_F182': 'linoleic_8', 'NUT_F160': 'palmitic_8'},
                       axis='columns')
               # the below could be written more elegantly
               .loc[:, lambda x: (x.columns == 'subjID') | 
                    x.columns.str.contains('_8')]
               .assign(sfa_pct_8 = lambda x: x.sfa_8 * 9 / x.tot_cal_8, 
                       mufa_pct_8 = lambda x: x.mufa_8 * 9 / x.tot_cal_8, 
                       pufa_pct_8 = lambda x: x.pufa_8 * 9 / x.tot_cal_8, 
                       sfa2pufa_8 = lambda x: x.sfa_8 / x.pufa_8,
                       palmitic_pct_8 = lambda x: x.palmitic_8 * 9 / 
                       x.tot_cal_8,
                       tot_fat_8 = lambda x: x.sfa_8 + x.mufa_8 + x.pufa_8,
                       tot_fat_pct_8 = lambda x: x.tot_fat_8 / x.tot_cal_8,
                       f2c_8 = lambda x: x.tot_fat_8 / x.carb_8,
                       n62n3_8 = lambda x: x.n6_8 / x.n3_8))

ffq_fhs = pd.merge(ffq_fhs_ex5, ffq_fhs_ex8, on="subjID", how="outer")

# Merge and save
               
fhs_metadata = pd.merge(phenos_fhs, ffq_fhs, on="subjID")
fhs_metadata.to_csv("../data/processed/metadata_fhs.csv", index=False)

# GWAS prep

# Pedigree
ped = (pd.read_csv("../data/raw/fhs/fhs_pedigree.txt", sep="\t", skiprows=10)
       .loc[:, ['pedno', 'shareid']]
       .rename(columns={'pedno': 'FID', 'shareid': 'IID'}))

# Construct outcome feature and prepare plink-friendly dataset
def rank_INT(values):
    # Computes rank-based inverse normal transformation of a Pandas series
    c = 3.0 / 8
    ranks = stats.rankdata(values)
    product_INT = stats.norm.ppf((ranks - c) / (len(ranks) - 2*c + 1))
    return product_INT

gwas_variables = ['subjID', 'ldl_5', 'sfa_pct_5', 'pufa_pct_5', 'tot_fat_pct_5',
                  'carb_5', 'palmitic_pct_5', 'f2c_5', 'age_5', 'tg_5', 'glu_5', 
                  'sbp_5', 'bmi_5', 'ldl_resid_5', 'sfa_pct_resid_5', 
                   'sfa2pufa_5', 'n62n3_5']
phenotypes = ['sfa_ldl_5', 'fat_ldl_5', 'f2c_ldl_5', 'palm_ldl_5',
              'sfa_bmi_5', 'f2c_bmi_5', 'sfa_sbp_5', 'f2c_sbp_5',
              'sfa_ldlR_5', 'sfaR_ldl_5', 'sfaR_ldlR_5',
              'f2c_tg_5', 'f2c_glu_5',
              'sfa2pufa_ldlR_5', 'n62n3_ldlR_5', 'n62n3_tg_5',
              'sfa2pufa_tg_5']

fhs_metadata_nomeds = fhs_metadata.query('lipid_med_5 == False').copy()
fhs_metadata_nomeds["ldl_resid_5"] = smf.ols('ldl_5 ~ age_5 + sex + bmi_5', 
                                             data=fhs_metadata_nomeds).fit().resid
fhs_metadata_nomeds["sfa_pct_resid_5"] = smf.ols('sfa_pct_5 ~ pufa_pct_5', 
                                                 data=fhs_metadata_nomeds).fit().resid

gwas_phenos = (fhs_metadata_nomeds
               .query('lipid_med_5 == False')
               .filter(gwas_variables + ['sex'])
               .merge(ped, left_on="subjID", right_on="IID") 
               .dropna()
               .assign(FID = lambda x: x.subjID,
                       IID = lambda x: x.subjID)
               .assign(sfa_ldl_5 = lambda x: scale(x.sfa_pct_5) * scale(x.ldl_5),
                       sfa_ldlR_5 = lambda x: scale(x.sfa_pct_5) *
                       scale(x.ldl_resid_5),
                       sfaR_ldl_5 = lambda x: scale(x.sfa_pct_resid_5) *
                       scale(x.ldl_5),
                       sfaR_ldlR_5 = lambda x: scale(x.sfa_pct_resid_5) *
                       scale(x.ldl_resid_5),
                       fat_ldl_5 = lambda x: scale(x.tot_fat_pct_5) * scale(x.ldl_5),
                       f2c_ldl_5 = lambda x: scale(x.f2c_5) * scale(x.ldl_5),
                       palm_ldl_5 = lambda x: scale(x.palmitic_pct_5) *
                       scale(x.ldl_5),
                       sfa_bmi_5 = lambda x: scale(x.sfa_pct_5) * scale(x.bmi_5),
                       f2c_bmi_5 = lambda x: scale(x.f2c_5) * scale(x.bmi_5),
                       f2c_tg_5 = lambda x: scale(x.f2c_5) * scale(x.tg_5),
                       f2c_glu_5 = lambda x: scale(x.f2c_5) * scale(x.glu_5),
                       sfa_sbp_5 = lambda x: scale(x.sfa_pct_5) *
                       scale(x.sbp_5),
                       f2c_sbp_5 = lambda x: scale(x.f2c_5) * scale(x.sbp_5),
                       sfa2pufa_5 = lambda x: x.sfa_pct_5 / x.pufa_pct_5,
                       sfa2pufa_ldlR_5 = lambda x: scale(x.sfa2pufa_5) *
                       scale(x.ldl_resid_5),
                       sfa2pufa_tg_5 = lambda x: scale(x.sfa2pufa_5) *
                       scale(x.tg_5),
                       n62n3_ldlR_5 = lambda x: scale(x.n62n3_5) *
                       scale(x.ldl_resid_5),
                       n62n3_tg_5 = lambda x: scale(x.n62n3_5) *
                       scale(x.tg_5))
               .filter(['FID', 'IID', 'sex'] + phenotypes + gwas_variables))
               ## FOR THE MOMENT, HAVE NOT FILTERED OUT RELATED SUBJECTS
gwas_phenos_WIN = (gwas_phenos
                   .filter(phenotypes)
                   .apply(stats.mstats.winsorize, limits=[0.01, 0.01])
                   .rename({p: p + "_WIN" for p in phenotypes}, axis=1))
gwas_phenos_INT = (gwas_phenos
                   .filter(phenotypes)
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
               
gwas_phenos.to_csv("../data/processed/fhs/fhs_gwas_phenos.txt", sep=" ", 
                   na_rep="NA", index=False)
