import sys
import pandas as pd
import numpy as np
import statsmodels.formula.api as smf
from sklearn.preprocessing import scale


diet_var = sys.argv[1]
pheno = sys.argv[2]
outdir = sys.argv[3]

meta_whi = pd.read_csv(f"{outdir}/meta_whi.txt", sep=" ")

scores = (pd.read_csv(f"{outdir}/interaction_scores.profile", delim_whitespace=True)
          .assign(IID=lambda x: pd.to_numeric(x.IID, errors="coerce"))
          .loc[:,['IID','SCORE']] .dropna()
          .assign(IID=lambda x: x.IID.astype(int)))

test_ids = (pd.read_csv(f"{outdir}/whi_test_ids.txt", 
                        sep=" ", header=None, usecols=[0], names=["IID"])
            .assign(IID=lambda x: pd.to_numeric(x.IID, errors="coerce"))
            .dropna()
            .assign(IID=lambda x: x.IID.astype(int)))

regData = (pd.merge(meta_whi, scores, on="IID")
           .assign(INT=lambda x: x[diet_var] * x.SCORE,
                   TERT=lambda x: pd.qcut(x.SCORE, 3, labels=[1,2,3]))
           .loc[lambda x: np.isin(x.IID, test_ids.IID),:])

with open(f"{outdir}/test_res.txt", "w") as f:
    for tert in range(1,4):
        subset = regData[regData.TERT == tert].copy()
        subset.loc[:,diet_var] = scale(subset.loc[:,diet_var])
        res = smf.ols(f"{pheno} ~ {diet_var}", data=subset).fit()
        f.write(f"TERTILE {tert}:\n")
        res.df = pd.DataFrame({'COEF': res.params,
                               'P-VALUE': res.pvalues})
        res.df.loc[diet_var].to_string(f)
        f.write("\n\n")
