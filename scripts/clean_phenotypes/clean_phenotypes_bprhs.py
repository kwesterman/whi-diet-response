import sys
import numpy as np
import pandas as pd
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
                     palmitic_pct = lambda x: x.palmitic * 9 / x.tot_cal))

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
gwas_features = ['subjID', 'ldl', 'sfa_pct', 'pufa_pct', 'age', 'bmi', 'PC1']
gwas_phenos = (bprhs_metadata
               .query('lipid_med == False')
               .filter(gwas_features + ['sex'])
               .dropna()
               .assign(sfa_ldl_product = lambda x: scale(x.sfa_pct) *
                       scale(x.ldl))
               .merge(pc_df, on="subjID")
               .assign(FID = lambda x: x.nelid_b,
                       IID = lambda x: x.nelid_b)
               .filter(['FID', 'IID', 'sex', 'sfa_ldl_product'] +
                       gwas_features))
               
gwas_phenos.to_csv("../data/processed/bprhs/bprhs_gwas_phenos.txt", sep=" ", 
                   index=False)
