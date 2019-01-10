import sys
import numpy as np
import pandas as pd
from scipy import stats
from sas7bdat import SAS7BDAT
from sklearn.preprocessing import scale


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

def winsorize(values, SDs_from_mean=4):
    # Winsorize values based on SDs from mean (rather than percentiles)
    sd = np.std(values)
    values[values < values.mean() - 4 * sd] = values.mean() - 4 * sd
    values[values > values.mean() + 4 * sd] = values.mean() + 4 * sd
    return values

gwas_variables = ['subjID', 'ldl', 'sfa_pct', 'pufa_pct', 'tot_fat_pct', 'carb',
                  'f2c', 'age', 'bmi', 'tg', 'glu', 'PC1']
phenotypes = ['sfa_ldl', 'fat_ldl', 'f2c_ldl', 
              'sfa_bmi', 'f2c_bmi', 'f2c_tg', 'f2c_glu']

gwas_phenos = (bprhs_metadata
               .query('lipid_med == False')
               .filter(gwas_features + ['sex'])
               .dropna()
               .merge(pc_df, on="subjID") 
               .assign(FID = lambda x: x.nelid_b, IID = lambda x: x.nelid_b) 
               .assign(sfa_ldl = lambda x: scale(x.sfa_pct) * scale(x.ldl),
                       fat_ldl = lambda x: scale(x.tot_fat_pct) * scale(x.ldl),
                       f2c_ldl = lambda x: scale(x.f2c) * scale(x.ldl),
                       sfa_bmi = lambda x: scale(x.sfa_pct) * scale(x.bmi),
                       f2c_bmi = lambda x: scale(x.f2c) * scale(x.bmi),
                       f2c_tg = lambda x: scale(x.f2c) * scale(x.tg),
                       f2c_glu = lambda x: scale(x.f2c) * scale(x.glu))
               .filter(['FID', 'IID', 'sex'] + phenotypes + gwas_variables))
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
               
gwas_phenos.to_csv("../data/processed/bprhs/bprhs_gwas_phenos.txt", sep=" ", 
                   index=False)
