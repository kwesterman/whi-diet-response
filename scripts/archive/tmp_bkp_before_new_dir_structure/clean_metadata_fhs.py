import sys
import pandas as pd
import statsmodels.formula.api as smf
from sklearn.preprocessing import scale


outdir = sys.argv[1]

phenos = (pd.read_csv("../int/metaData_fhs.csv")
          .query('lipid_med_5 == False')
          .loc[:, ['subjID', 'ldl_5', 'sfa_5', 'pufa_5', 'age_5', 'sex']]
          .dropna())

phenos['ldl_adj_5'] = smf.ols('ldl_5 ~ pufa_5 + age_5 + sex', data=phenos).fit().resid
phenos['int_product'] = scale(phenos.sfa_5) * scale(phenos.ldl_adj_5)

first_cols = ['subjID', 'subjID', 'sex']
phenos = phenos.loc[:, lambda df: first_cols + list(df.columns.difference(first_cols))]
phenos.columns = ['FID', 'IID'] + list(phenos.columns[2:])

phenos.to_csv(f"{outdir}/meta_fhs.txt", sep=" ", index=False)
phenos[['FID', 'IID']].to_csv(f"{outdir}/train_ids_fhs.txt", sep=" ",
                              index=False)
