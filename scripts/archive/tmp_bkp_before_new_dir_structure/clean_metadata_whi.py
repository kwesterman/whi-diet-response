import sys
import pandas as pd
import statsmodels.formula.api as smf
from sklearn.preprocessing import scale


outdir = sys.argv[1]

# Sample-to-subject mapping
sample_ids = (pd.read_csv("../data/whi/gen/imputedSampleInfo.txt", sep="\t", skiprows=15)
              .loc[:,['SampleID','SubjectID']])

# Contains ancestry principal components linked to samples
sample_info_c1 = (pd.read_csv("../data/whi/gen/sampleInfoWithPCs_c1.txt", sep="\t", skiprows=10)
                  .loc[:,['SAMPLE_ID','PC1','PC2','PC3','PC4','PC5']])
sample_info_c2 = (pd.read_csv("../data/whi/gen/sampleInfoWithPCs_c2.txt", sep="\t", skiprows=10)
                  .loc[:,['SAMPLE_ID','PC1','PC2','PC3','PC4','PC5']])
sample_info = pd.concat([sample_info_c1, sample_info_c2], axis=0).reset_index()

pc_df = (pd.merge(sample_ids, sample_info, left_on='SampleID', right_on='SAMPLE_ID')
         .filter(regex="SubjectID|PC.*")
         .groupby('SubjectID')
         .first())

# Main diet and phenotype dataset
phenos = (pd.read_csv("../int/metaData_whi.csv")
          .query('lipid_med == False and visitYear == 0 and dm_trial == False')
          .loc[:, ['subjID', 'ldl', 'sfa', 'pufa', 'age', 'sex']]
          .dropna())

phenos['ldl_adj'] = smf.ols('ldl ~ pufa + age', data=phenos).fit().resid
phenos['int_product'] = scale(phenos.sfa) * scale(phenos.ldl_adj)

phenos = phenos.merge(pc_df, left_on="subjID", right_index=True)

first_cols = ['subjID', 'subjID', 'sex']
phenos = phenos[first_cols + list(phenos.columns.difference(first_cols))]
phenos.columns = ['FID','IID'] + list(phenos.columns)[2:]

phenos.to_csv(f"{outdir}/meta_whi.txt", sep=" ", index=False)
phenos[['FID', 'IID']].to_csv(f"{outdir}/train_ids_whi.txt", sep=" ",
                              index=False)
