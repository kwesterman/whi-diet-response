import sys
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.preprocessing import scale
import statsmodels.formula.api as smf
import itertools


whi_metadata = pd.read_csv("../../data/processed/metadata_whi.csv")

# Sample-to-subject mapping
sample_ids = (pd.read_csv("../../data/raw/whi/gen/imputedSampleInfo.txt", sep="\t", skiprows=15)
              .rename({'SubjectID': 'subjID'}, axis=1)
              .filter(['SampleID', 'subjID']))

# Contains ancestry principal components linked to samples
sample_info_c1 = (pd.read_csv("../../data/raw/whi/gen/sampleInfoWithPCs_c1.txt", sep="\t", skiprows=10)
                  .filter(['SAMPLE_ID','PC1','PC2','PC3','PC4','PC5']))
sample_info_c2 = (pd.read_csv("../../data/raw/whi/gen/sampleInfoWithPCs_c2.txt", sep="\t", skiprows=10)
                  .filter(['SAMPLE_ID','PC1','PC2','PC3','PC4','PC5']))
sample_info = pd.concat([sample_info_c1, sample_info_c2], ignore_index=True)

pc_df = (pd.merge(sample_ids, sample_info, left_on='SampleID', right_on='SAMPLE_ID')
         .filter(regex="subjID|SampleID|PC.*")
         .groupby('subjID')
         .first()  # Arbitrarily, takes the "first" sample for each subject
         .reset_index())

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

def winsorize(x, num_sd=5):
    return x.clip(lower=x.mean() - num_sd * x.std(),
                  upper=x.mean() + num_sd * x.std())


gwas_covars = ['visit_year', 'sex', 'race', 'age', 'lipid_med', 'ht_med', 'dm_med'] + ["PC" + str(i) for i in range(1, 6)]

diet_vars = ['sfa', 'mufa', 'pufa', 'fat', 'palm', 'n3', 'pro', 'cho', 'tot_cal', 'alc', 'myhei', 'hei']
             
risk_factors = ["chol", "ldl", "hdl", "tg", "glu", "sbp", "bmi", 'hsCRP', 'logHDL', 'logTG', 'logBMI', 'logGLU', 'logSBP', 'logHSCRP']

deltas = (whi_metadata
	  .query('dm_trial == False')
          .query('subjID in @pc_df.subjID')
          .filter(['subjID', 'visit_year'] + diet_vars + risk_factors) 
          .groupby('subjID')
          .filter(lambda x: x.shape[0] > 1)
          .sort_values('visit_year')
          .groupby('subjID')
          .agg(lambda x: x.iloc[-1] - x.iloc[0])
          .reset_index()
          .drop('visit_year', axis=1)
          .rename({var: "delta_" + var for var in diet_vars + risk_factors},
                  axis=1)).copy()

binaries = (whi_metadata
	    .query('dm_trial == False')
	    .query('visit_year == 0')
	    .query('subjID in @pc_df.subjID')
	    .filter(['subjID'] + diet_vars))
binaries.loc[:, diet_vars] = (binaries.loc[:, diet_vars]
			      #.apply(lambda x: pd.qcut(x, [0, .25, .5, .75, 1.], labels=False))
			      .apply(lambda x: pd.qcut(x, [0, 0.333, 0.666, 1.], labels=False))
			      .apply(lambda x: np.where(x.isin([0., 1.]), x, np.nan)))
binaries = binaries.rename({dv: dv + "Bin" for dv in diet_vars}, axis=1)

diet_vars = diet_vars + ['sfa_binary', 'sfa_pct_binary', 'fat_binary', 'fat_pct_binary']

whi_metadata = (whi_metadata
		.assign(ldl = lambda x: np.where(x.lipid_med, x.ldl / 0.75, x.ldl),
			glu = lambda x: np.where(x.dm_med, x.glu / 0.75, x.glu),
			sbp = lambda x: np.where(x.sbp, x.sbp + 15, x.sbp)))

gwas_phenos = (whi_metadata
               .query('dm_trial == False')
               .query('visit_year == 0')
               .filter(['subjID'] + gwas_covars + diet_vars + risk_factors)
               #.dropna()
               .drop_duplicates(subset=["subjID"])
               .merge(pc_df, on="subjID")
	       .merge(deltas, on="subjID", how="left")
	       .merge(binaries, on="subjID", how="left"))

diet_vars.extend(["delta_" + dv for dv in diet_vars])
risk_factors.extend(["delta_" + rf for rf in risk_factors])
diet_vars.extend([dv + "Bin" for dv in diet_vars])

gwas_phenos_noFilter = (gwas_phenos
	       		.assign(ldl=lambda x: winsorize(x.ldl),
				logHDL=lambda x: np.log(winsorize(x.hdl)),
			        logTG=lambda x: np.log(winsorize(x.tg)),
			        logBMI=lambda x: np.log(winsorize(x.bmi)),
			        logGLU=lambda x: np.log(winsorize(x.glu)),
			        logHSCRP=lambda x: np.log(winsorize(x.hsCRP)),
			        logSBP=lambda x: np.log(winsorize(x.sbp)),
				sfa_binary=lambda x: (x.sfa > 22) * 1,
				sfa_pct=lambda x: x.sfa * 9 / x.tot_cal * 100,
				sfa_pct_binary=lambda x: (x.sfa_pct > x.sfa_pct.median()) * 1,
				fat_binary=lambda x: (x.fat > 22) * 1,
				fat_pct=lambda x: x.fat * 9 / x.tot_cal * 100,
				fat_pct_binary=lambda x: (x.fat_pct > x.fat_pct.median()) * 1)
               .assign(FID = lambda x: x.SampleID, IID = lambda x: x.SampleID)
               .filter(['FID', 'IID', 'sex'] + gwas_covars + diet_vars + risk_factors))

gwas_phenos_noFilter["ldl_resid"] = smf.ols('ldl ~ age + tot_cal', 
                                             data=gwas_phenos_noFilter).fit().resid
gwas_phenos_noFilter["hdl_resid"] = smf.ols('logHDL ~ age + tot_cal', 
                                             data=gwas_phenos_noFilter).fit().resid
gwas_phenos_noFilter["fat_ldl_prod"] = ((pd.to_numeric(gwas_phenos_noFilter.fat_pct_binary, errors="coerce") * 2) - 1) * gwas_phenos_noFilter.ldl_resid
gwas_phenos_noFilter["fat_hdl_prod"] = ((pd.to_numeric(gwas_phenos_noFilter.fat_pct_binary, errors="coerce") * 2) - 1) * gwas_phenos_noFilter.hdl_resid

gwas_phenos_noFilter.query('race == "white"').to_csv(
    "../../data/processed/gen6/whi_white_gwas_phenos.txt", sep=" ", na_rep="NA", 
    index=False)
#gwas_phenos.query('race == "black"').to_csv(
#    "../../data/processed/gen6/whi_black_gwas_phenos_long.txt", sep=" ", na_rep="NA", 
#    index=False)
#gwas_phenos.query('race == "hispanic"').to_csv(
#    "../../data/processed/gen6/whi_hispanic_gwas_phenos_long.txt", sep=" ", na_rep="NA", 
#    index=False)


#gwas_phenos_noLipidMed = (gwas_phenos
#	       .query('lipid_med == False')
#	       .assign(ldl=lambda x: winsorize(x.ldl),
#		       logTG=lambda x: np.log(winsorize(x.tg)),
#		       logBMI=lambda x: np.log(winsorize(x.bmi)),
#		       logGLU=lambda x: np.log(winsorize(x.glu)),
#		       logHSCRP=lambda x: np.log(winsorize(x.hsCRP)),
#		       logSBP=lambda x: np.log(winsorize(x.sbp)))
#               .assign(FID = lambda x: x.SampleID, IID = lambda x: x.SampleID)
#               .filter(['FID', 'IID', 'sex'] + gwas_covars + diet_vars + risk_factors))
#gwas_phenos_noBPMed = (gwas_phenos
#	       .query('ht_med == False')
#	       .assign(ldl=lambda x: winsorize(x.ldl),
#		       logTG=lambda x: np.log(winsorize(x.tg)),
#		       logBMI=lambda x: np.log(winsorize(x.bmi)),
#		       logGLU=lambda x: np.log(winsorize(x.glu)),
#		       logHSCRP=lambda x: np.log(winsorize(x.hsCRP)),
#		       logSBP=lambda x: np.log(winsorize(x.sbp)))
#               .assign(FID = lambda x: x.SampleID, IID = lambda x: x.SampleID)
#               .filter(['FID', 'IID', 'sex'] + gwas_covars + diet_vars + risk_factors))
#gwas_phenos_noDMMed = (gwas_phenos
#	       .query('dm_med == False')
#	       .assign(ldl=lambda x: winsorize(x.ldl),
#		       logTG=lambda x: np.log(winsorize(x.tg)),
#		       logBMI=lambda x: np.log(winsorize(x.bmi)),
#		       logGLU=lambda x: np.log(winsorize(x.glu)),
#		       logHSCRP=lambda x: np.log(winsorize(x.hsCRP)),
#		       logSBP=lambda x: np.log(winsorize(x.sbp)))
#               .assign(FID = lambda x: x.SampleID, IID = lambda x: x.SampleID)
#               .filter(['FID', 'IID', 'sex'] + gwas_covars + diet_vars + risk_factors))
