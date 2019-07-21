import sys
import numpy as np
import pandas as pd
from scipy import stats
from sas7bdat import SAS7BDAT
from sklearn.preprocessing import scale
import statsmodels.formula.api as smf


# Phenotypes

phenos_list = pd.read_sas("../data/raw/bprhs/phen/csbl_26sep2013_list.sas7bdat")
phenos = (SAS7BDAT("../data/raw/bprhs/phen/csbl_26sep2013.sas7bdat")
          .to_data_frame()
          .assign(sex = lambda x: np.where(x.female == 1, "F", "M"))
          .rename({'BMI': 'bmi', 'trig': 'tg', 'gluc': 'glu', 'crp': 'hsCRP',
                   'sysbp': 'sbp', 'mhmg': 'lipid_med'}, axis=1)
          .filter(['studyid', 'sex', 'age', 'bmi', 'ldl', 'vldl', 'hdl', 'tg',
                   'glu', 'hsCRP', 'sbp', 'lipid_med']))

# FFQ

ffq_bprhs = (SAS7BDAT("../data/raw/bprhs/diet/cphhdv1b_datot_111512.sas7bdat")
             .to_data_frame()
             .rename({'kcal_f': 'tot_cal', 'sfa_f': 'sfa', 'mfa_f': 'mufa',
                      'pfa_f': 'pufa', 'tcho_f': 'carb', 's16_0_f': 'palmitic',
                      'p18_2_f': 'linoleic', 'omega3_fs': 'n3'}, axis=1)
             .filter(['studyid', 'tot_cal', 'sfa', 'mufa', 'pufa', 'carb',
                      'palmitic', 'linoleic', 'n3', 'n6'])
             .assign(n6 = lambda x: x.linoleic,  # Workaround for no explicit n6
                     sfa_pct = lambda x: x.sfa * 9 / x.tot_cal,
                     mufa_pct = lambda x: x.mufa * 9 / x.tot_cal,
                     pufa_pct = lambda x: x.pufa * 9 / x.tot_cal,
                     palmitic_pct = lambda x: x.palmitic * 9 / x.tot_cal,
                     tot_fat = lambda x: x.sfa + x.mufa + x.pufa,
                     tot_fat_pct = lambda x: x.tot_fat / x.tot_cal,
                     f2c = lambda x: x.tot_fat / x.carb))

bprhs_metadata = (pd.merge(phenos, ffq_bprhs, on="studyid")
                  .rename({'studyid': 'subjID'}, axis=1))
bprhs_metadata.to_csv("../data/processed/metadata_bprhs.csv", index=False)


# GWAS prep

# Sample-to-subject mapping
sample_ids = pd.read_csv("../data/raw/bprhs/admin/BPRHS_id_link.csv",
                         dtype={"studyid": "object", "nelid_b": "object"})

# Ancestry principal components linked to samples
pc_df = (pd.read_csv("../data/raw/bprhs/gen/bprhs_pca.csv", 
                     dtype={"nelid_b": "object"})
         .merge(sample_ids, on="nelid_b")
         .filter(['studyid', 'nelid_b', 'PCA1'])
         .rename({'studyid': 'subjID', 'PCA1': 'PC1'}, axis=1))
### FUTURE: THE PC DATA FRAME CONTAINS ONLY ABOUT HALF OF THE INDIVIDUALS

# Construct outcome "interaction" feature and prepare plink-friendly dataset
def rank_INT(values):
    # Computes rank-based inverse normal transformation of a Pandas series
    c = 3.0 / 8
    ranks = stats.rankdata(values)
    product_INT = stats.norm.ppf((ranks - c) / (len(ranks) - 2*c + 1))
    return product_INT

gwas_variables = ['subjID', 'ldl', 'sfa_pct', 'pufa_pct', 'tot_fat_pct', 'carb',
                  'palmitic_pct', 'f2c', 'age', 'bmi', 'tg', 'glu', 'sbp',
                  'ldl_resid', 'sfa_pct_resid', 'PC1']
phenotypes = ['sfa_ldl', 'fat_ldl', 'f2c_ldl', 'palm_ldl', 
              'sfa_ldlR', 'sfaR_ldl', 'sfaR_ldlR',
              'sfa_bmi', 'f2c_bmi', 'sfa_sbp', 'f2c_sbp', 'f2c_tg', 'f2c_glu']

bprhs_metadata_nomeds = bprhs_metadata.query('lipid_med == False').copy()
bprhs_metadata_nomeds["ldl_resid"] = smf.ols('ldl ~ age + sex + bmi', 
                                             data=bprhs_metadata_nomeds).fit().resid
bprhs_metadata_nomeds["sfa_pct_resid"] = smf.ols('sfa_pct ~ pufa_pct', 
                                                 data=bprhs_metadata_nomeds).fit().resid

gwas_phenos = (bprhs_metadata_nomeds
               .query('lipid_med == False')
               .filter(gwas_variables + ['sex'])
               .dropna()
               .merge(pc_df, on="subjID") 
               .assign(FID = lambda x: x.nelid_b, IID = lambda x: x.nelid_b) 
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
                       f2c_sbp = lambda x: scale(x.f2c) * scale(x.sbp))
               .filter(['FID', 'IID', 'sex'] + phenotypes + gwas_variables))
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
                       #~x.between(*x.quantile([0.25, 0.75]))))
                   .rename({(p + "_WIN"): p + "_BIN" for p in phenotypes}, 
                           axis=1))
gwas_phenos = pd.concat([gwas_phenos, gwas_phenos_WIN, gwas_phenos_INT,
                         gwas_phenos_BIN], axis=1) 
               
gwas_phenos.to_csv("../data/processed/bprhs/bprhs_gwas_phenos.txt", sep=" ", 
                   na_rep="NA", index=False)
