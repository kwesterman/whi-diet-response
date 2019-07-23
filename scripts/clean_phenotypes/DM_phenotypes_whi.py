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


diet_vars = ['sfa', 'mufa', 'pufa', 'fat', 'palm', 'n3', 'pro', 'cho', 'tot_cal', 'alc']
             
risk_factors = ["chol", "ldl", "hdl", "tg", "glu", "sbp", "bmi", 'hsCRP', 'logTG', 'logBMI', 'logGLU', 'logSBP', 'logHSCRP']

addtl_vars = ["sex", "race", "age", "lipid_med", "ht_med", "dm_med", "dm_trial", "dm_intervention"] + ["PC" + str(i) for i in range(1, 6)]

deltas = (whi_metadata
	  .query('dm_trial == True')
          .query('subjID in @pc_df.subjID')
          .filter(['subjID', 'visit_year'] + diet_vars + risk_factors) 
          #.dropna(subset=diet_vars + risk_factors)
          .groupby('subjID')
          .filter(lambda x: x.shape[0] > 1)
          .sort_values('visit_year')
          .groupby('subjID')
          .agg(lambda x: x.iloc[-1] - x.iloc[0])
          .reset_index()
          .drop('visit_year', axis=1)
          .rename({var: "delta_" + var for var in diet_vars + risk_factors},
                  axis=1))

baselines = (whi_metadata
	     .query('visit_year == 0')
	     .filter(["subjID"] + diet_vars + risk_factors))
#baselines.loc[:, baselines.columns != "subjID"] = baselines.loc[:, baselines.columns != "subjID"].add_prefix("baseline_")
baselines.columns = ['baseline_' + nm if nm != "subjID" else nm for nm in baselines.columns]

dm_long_phenos = (whi_metadata
		  .filter(["subjID"] + addtl_vars)
		  .drop_duplicates(subset="subjID")
		  .merge(pc_df, on="subjID")
		  .merge(baselines, how="left", on="subjID")
		  .merge(deltas, how="left", on="subjID"))

dm_long_phenos.query('race == "white"').to_csv(
    "../../data/processed/gen6/whi_white_DM_phenos.txt", sep=" ", na_rep="NA", 
    index=False)
dm_long_phenos.query('race == "black"').to_csv(
    "../../data/processed/gen6/whi_black_DM_phenos.txt", sep=" ", na_rep="NA", 
    index=False)
dm_long_phenos.query('race == "hispanic"').to_csv(
    "../../data/processed/gen6/whi_hispanic_DM_phenos.txt", sep=" ", na_rep="NA", 
    index=False)
