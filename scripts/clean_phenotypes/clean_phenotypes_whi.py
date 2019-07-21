import sys
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.preprocessing import scale
import statsmodels.formula.api as smf
import itertools


# Risk factors

basic_data_whi_c1 = pd.read_csv("../data/raw/whi/phen/basic_whi_c1.txt",
                                sep="\t", skiprows=10)
basic_data_whi_c2 = pd.read_csv("../data/raw/whi/phen/basic_whi_c2.txt",
                                sep="\t", skiprows=10)
basic_data_whi = (pd.concat([basic_data_whi_c1, basic_data_whi_c2],
                            ignore_index=True)
                  .rename({'SUBJID': 'subjID', 'AGE': 'age', 'RACE': 'race'},
                          axis=1)
                  .assign(sex = "F")
                  .replace({1: 'amind', 2: 'asian', 3: 'black', 4: 'hispanic',
                            5: 'white', 8: 'other'})
                  .filter(['subjID', 'age', 'sex', 'race']))


randomization_whi_c1 = pd.read_csv("../data/raw/whi/randomization_c1.txt",
                                sep="\t", skiprows=10)
randomization_whi_c2 = pd.read_csv("../data/raw/whi/randomization_c2.txt",
                                sep="\t", skiprows=10)
randomization_whi = (pd.concat([randomization_whi_c1, randomization_whi_c2],
                               ignore_index=True)
                     .rename({'SUBJID': 'subjID'}, axis=1)
                     .assign(sex = "F")
                     .assign(dm_trial = lambda x: x.DMFLAG.astype('bool'),
                             dm_intervention = lambda x: x.DMARM == 1)
                     .filter(['subjID', 'dm_trial', 'dm_intervention']))


behavior_data_whi_c1 = pd.read_csv("../data/raw/whi/phen/behavior_c1.txt",
                                sep="\t", skiprows=10)
behavior_data_whi_c2 = pd.read_csv("../data/raw/whi/phen/behavior_c2.txt",
                                sep="\t", skiprows=10)
behavior_data_whi = (pd.concat([behavior_data_whi_c1, behavior_data_whi_c2],
                               ignore_index=True)
                     .rename({'SUBJID': 'subjID', 'SMOKNOW': 'smk_now',
                              'PACKRYS': 'smk_py'}, axis=1)
                     .filter(['subjID', 'smk_now', 'smk_py']))

exam_data_whi_c1 = pd.read_csv("../data/raw/whi/phen/physical_whi_c1.txt",
                                sep="\t", skiprows=10)
exam_data_whi_c2 = pd.read_csv("../data/raw/whi/phen/physical_whi_c2.txt",
                                sep="\t", skiprows=10)
exam_data_whi = (pd.concat([exam_data_whi_c1, exam_data_whi_c2],
                           ignore_index=True)
                     .rename({'SUBJID': 'subjID', 'F80VY': 'visit_year',
                              'SYST': 'sbp', 'BMIX': 'bmi'}, axis=1)
                     .filter(['subjID', 'visit_year', 'sbp', 'bmi']))


lab_data1_whi_c1 = pd.read_csv("../data/raw/whi/phen/labs_c1.txt",
                                sep="\t", skiprows=10)
lab_data1_whi_c2 = pd.read_csv("../data/raw/whi/phen/labs_c2.txt",
                                sep="\t", skiprows=10)
lab_data1_whi = (pd.concat([lab_data1_whi_c1, lab_data1_whi_c2],
                           ignore_index=True)
                 .rename({'SUBJID': 'subjID', 'COREVY': 'visit_year',
                          'COREGLUC': 'glu', 'CORETCHO': 'chol',
                          'CORELDLC': 'ldl', 'COREHDLC': 'hdl', 
                          'CORETRI': 'tg'}, axis=1)
                 .filter(['subjID', 'visit_year', 'glu', 'chol', 'ldl',
                          'hdl', 'tg']))


draw_data2_whi_c1 = pd.read_csv("../data/raw/whi/phen/draws2_c1.txt", 
                               sep="\t", skiprows=10, 
                               usecols=['DRAWID', 'DRAWVTYP', 'DRAWVY'])
draw_data2_whi_c2 = pd.read_csv("../data/raw/whi/phen/draws2_c2.txt",
                               sep="\t", skiprows=10,
                               usecols=['DRAWID', 'DRAWVTYP', 'DRAWVY'])
draw_data2_whi = (pd.concat([draw_data2_whi_c1, draw_data2_whi_c2],
                            ignore_index=True)
                  .query('DRAWVTYP != 4')  # Throw out non-routine visits
                  .rename({'DRAWVY': 'visit_year'}, axis=1)
                  .filter(['DRAWID', 'visit_year']))

lab_data2_whi_c1 = pd.read_csv("../data/raw/whi/phen/labs2_c1.txt", 
                               sep="\t", skiprows=10, encoding = "ISO-8859-1")
lab_data2_whi_c2 = pd.read_csv("../data/raw/whi/phen/labs2_c2.txt",
                               sep="\t", skiprows=10, encoding = "ISO-8859-1")
lab_data2_whi = (pd.concat([lab_data2_whi_c1, lab_data2_whi_c2],
                           ignore_index=True)
                 .query('SPECTYPE == "Serum"')
                 .loc[lambda x: x.TESTABBR.isin(['LDLC', 'HDLC', 'TCHO', 'TRI',
                                                'GLUC', 'CRP', 'INSU'])]
                 .merge(draw_data2_whi, on="DRAWID")
                 .groupby(['SUBJID', 'visit_year', 'TESTABBR'])
                 .agg({'TESTVAL': 'mean'})
                 .unstack(level='TESTABBR'))
lab_data2_whi.columns = lab_data2_whi.columns.droplevel(0)
lab_data2_whi.columns.name = None
lab_data2_whi = (lab_data2_whi
                 .reset_index()
                 .rename({'SUBJID': 'subjID', 'TCHO': 'chol', 'LDLC': 'ldl',
                          'HDLC': 'hdl', 'TRI': 'tg', 'GLUC': 'glu',
                          'CRP': 'hsCRP', 'INSU': 'ins'}, axis=1))

lab_data_whi = (pd.concat([lab_data1_whi, lab_data2_whi], sort=False,
                          ignore_index=True)
                .groupby(['subjID', 'visit_year'])
                .agg('mean'))

meds_data_whi_c1 = pd.read_csv("../data/raw/whi/phen/medications_c1.txt",
                               sep="\t", skiprows=10)
meds_data_whi_c2 = pd.read_csv("../data/raw/whi/phen/medications_c2.txt",
                               sep="\t", skiprows=10)
meds_ref_whi = pd.read_csv("../data/raw/whi/phen/medication_classes.dat",
                           sep="\t") 
lipid_med_classes = [
    'ANTIHYPERLIPIDEMIC', 'FIBRIC ACID DERIVATIVES', 
    'INTESTINAL CHOLESTEROL ABSORPTION INHIBITORS', 'HMG COA REDUCTASE INHIBITORS', 
    'HMG COA REDUCTASE INHIBITOR COMBINATIONS', 'NICOTINIC ACID DERIVATIVES', 
    'MISC. ANTIHYPERLIPIDEMICS', 'ANTIHYPERLIPIDEMIC COMBINATIONS', 
    'FIBRIC ACID DERIVATIVE COMBINATIONS', 
    'CALCIUM BLOCKER & HMG COA REDUCTASE INHIBITOR COMB']
diabetes_med_classes = [
    "ANTIDIABETIC", "INSULIN", "MIXED INSULIN", "BEEF INSULIN", "PORK INSULIN", 
    "HUMAN INSULIN", "SULFONYLUREAS", "SULFONYLUREA COMBINATIONS", 
    "ANTIDIABETIC - AMINO ACID DERIVATIVES", "ANTIDIABETIC - D-PHENYLALANINE DERIVATIVES", 
    "BIGUANIDES", "MEGLITINIDE ANALOGUES", "DIABETIC OTHER", 
    "DIABETIC OTHER - COMBINATIONS", "ALPHA-GLUCOSIDASE INHIBITORS", 
    "INSULIN SENSITIZING AGENTS", "THIAZOLIDINEDIONES", "ANTIDIABETIC COMBINATIONS", 
    "SULFONYLUREA-BIGUANIDE COMBINATIONS", 
    "THIAZOLIDINEDIONE-BIGUANIDE COMBINATIONS"]
hypertension_med_classes = [
    "MINERALOCORTICOIDS", "VASOPRESSIN", "BETA BLOCKERS", "BETA BLOCKERS NON-SELECTIVE", 
    "BETA BLOCKERS CARDIO-SELECTIVE", "ALPHA-BETA BLOCKERS", "CALCIUM BLOCKERS", 
    "ANTIHYPERTENSIVE", "ACE INHIBITORS", "ANGIOTENSIN II RECEPTOR ANTAGONIST", 
    "ANTIADRENERGIC ANTIHYPERTENSIVES", "ANTIADRENERGICS - CENTRALLY ACTING", 
    "ANTIADRENERGICS - PERIPHERALLY ACTING", "RESERPINE", 
    "SELECTIVE ALDOSTERONE RECEPTOR ANTAGONISTS (SARAS)", "VASODILATORS", 
    "FLUOROQUINOLONE VASODILATORS", "DOPAMINE D1 RECEPTOR AGONISTS", "ANTIHYPERTENSIVE - MAOIS", 
    "MISC. ANTIHYPERTENSIVES", "ANTIHYPERTENSIVE COMBINATIONS", "RESERPINE COMBINATIONS", 
    "ACE INHIBITORS & CALCIUM BLOCKERS", "ACE INHIBITORS & THIAZIDE/THIAZIDE-LIKE", 
    "BETA BLOCKER & DIURETIC COMBINATIONS", "BETA BLOCKER & CALCIUM BLOCKER COMBINATIONS", 
    "ANGIOTENSIN II RECEPTOR ANTAGONISTS & THIAZIDES", 
    "ADRENOLYTICS-CENTRAL & THIAZIDE COMBINATIONS", "ADRENOLYTICS-PERIPHERAL & THIAZIDES", 
    "ANTIHYPERTENSIVES-MAOIS & THIAZIDES", "ANTIHYPERTENSIVES-MISC & THIAZIDES", 
    "VASODILATORS & THIAZIDES", "DIURETICS", "CARBONIC ANHYDRASE INHIBITORS", "LOOP DIURETICS", 
    "MERCURIAL DIURETICS", "OSMOTIC DIURETICS", "POTASSIUM SPARING DIURETICS", 
    "THIAZIDES AND THIAZIDE-LIKE DIURETICS", "MISCELLANEOUS DIURETICS", "COMBINATION DIURETICS", 
    "DIURETICS & POTASSIUM", "NON PRESCRIPTION DIURETICS", "PRESSORS", "PRESSOR COMBINATIONS", 
    "PERIPHERAL VASODILATORS", "VASODILATOR COMBINATIONS", "MICROVASODILATORS", 
    "PROSTAGLANDIN VASODILATORS", "VASOACTIVE NATRIURETIC PEPTIDES", "VASOCONSTRICTOR INHIBITORS", 
    "CALCIUM BLOCKER & HMG COA REDUCTASE INHIBITOR COMB"]
meds_data_whi = (pd.concat([meds_data_whi_c1, meds_data_whi_c2],
                           ignore_index=True)
                 #.query('F44VY == 1')
                 .merge(meds_ref_whi, on="TCCODE")
                 .assign(ht_med = lambda x:
                         x.TCNAME.isin(hypertension_med_classes),
                         lipid_med = lambda x:
                         x.TCNAME.isin(lipid_med_classes),
                         dm_med = lambda x:
                         x.TCNAME.isin(diabetes_med_classes))
                 .rename({'SUBJID': 'subjID', 'F44VY': 'visit_year'}, axis=1)
                 .groupby(['subjID', 'visit_year'])
                 .agg({'ht_med': np.any, 'lipid_med': np.any, 'dm_med': np.any})
                 .reset_index()
                 .filter(['subjID', 'visit_year', 'ht_med', 'lipid_med',
                          'dm_med']))
retro_meds_data_whi = (meds_data_whi
                       .query('visit_year == 1')
                       .assign(visit_year = 0))  # Assumes that meds info from
                                                 # visit 1 also applied at baseline
meds_data_whi = (pd.concat([meds_data_whi, retro_meds_data_whi],
                           ignore_index=True)
                 .drop_duplicates())

misc_merged = (exam_data_whi
               .merge(lab_data_whi, on=["subjID", "visit_year"], how="outer")
               .merge(meds_data_whi, on=["subjID", "visit_year"], how="outer"))

pheno_data_whi = (basic_data_whi
                  .merge(randomization_whi, on="subjID", how="left")
                  .merge(behavior_data_whi, on="subjID", how="left")
                  .merge(misc_merged, on="subjID", how="left"))

# Diet

def winsorize(x, num_sd=5):
    return x.clip(lower=x.mean() - num_sd * x.std(),
                  upper=x.mean() + num_sd * x.std())

ffq_nutrients_whi_c1 = pd.read_csv("../data/raw/whi/diet/ffq_nutrients_c1.txt",
                                   sep="\t", skiprows=10)
ffq_nutrients_whi_c2 = pd.read_csv("../data/raw/whi/diet/ffq_nutrients_c2.txt",
                                   sep="\t", skiprows=10)
ffq_nutrients_whi = (pd.concat([ffq_nutrients_whi_c1, ffq_nutrients_whi_c2],
                               ignore_index=True)
                     .rename({'SUBJID': 'subjID', 'F60VY': 'visit_year', 
                              'F60ENRGY': 'tot_cal', 'F60FAT': 'fat', 
                              'F60PROT': 'pro', 'F60CARB': 'cho',
                              'F60SFA': 'sfa', 'F60MFA': 'mufa', 'F60PFA': 'pufa', 
                              'F60ALC': 'alc', 'F60FIBER': 'fiber',
                              'F60FRUIT': 'fruit', 'F60VEG': 'veg',
                              'F60SF160': 'palm', 'F60PF182': 'linoleic', 
                              'F60OMGA3': 'n3', 'F60OMGA6': 'n6', 'FRUVEG':
                              'FV', 'F60SODUM': 'sodium'}, axis=1)
                     .filter(['subjID', 'visit_year', 'tot_cal', 'fat', 'pro',
                              'cho', 'sfa', 'mufa', 'pufa', 'alc', 'fiber',
                              'fruit', 'veg', 'palm', 'linoleic', 'n3', 'n6', 
                              'FV', 'sodium'])
                     .assign(sfa_pct = lambda x: x.sfa * 9 / x.tot_cal,
                             mufa_pct = lambda x: x.mufa * 9 / x.tot_cal,
                             pufa_pct = lambda x: x.pufa * 9 / x.tot_cal,
                             logsfa2pufa = lambda x: np.log(x.sfa / x.pufa),
                             palm_pct = lambda x: x.palm * 9 / x.tot_cal,
                             fat_pct = lambda x: x.fat / x.tot_cal,
                             cho_pct = lambda x: x.cho / x.tot_cal,
                             logn62n3 = lambda x: np.log(x.n6 / x.n3),
                             logf2c = lambda x: np.log(x.fat / x.cho)))

ffq_items_whi_c1 = pd.read_csv("../data/raw/whi/diet/ffq_items_c1.txt",
                                   sep="\t", skiprows=10)
ffq_items_whi_c2 = pd.read_csv("../data/raw/whi/diet/ffq_items_c2.txt",
                                   sep="\t", skiprows=10)
ffq_items_whi = (pd.concat([ffq_items_whi_c1, ffq_items_whi_c2],
                               ignore_index=True)
                 .rename({'SUBJID': 'subjID', 'F60VY': 'visit_year', 
                          'REDMEAT': 'RM', 'WHLGRNS': 'WG', 'NUTS': 'nuts', 
                          'FRUITS': 'fruit', 'VEGTABLS': 'veg', 'FISH': 'fish',
                          'DAIRY': 'dairy', 'FRUVEG': 'FV', 'GRAINS': 'grains'}, axis=1)
                 .filter(['subjID', 'visit_year', 'FV', 'RM', 'WG', 'nuts',
                          'grains', 'dairy', 'fish'])
                 .assign(RM=lambda x: scale(winsorize(x.RM)),
                         WG=lambda x: scale(winsorize(x.WG)),
                         nuts=lambda x: scale(winsorize(x.nuts)),
                         FV=lambda x: scale(winsorize(x.FV)),
                         fish=lambda x: scale(winsorize(x.fish)),
                         dairy=lambda x: scale(winsorize(x.dairy)),
                         grains=lambda x: scale(winsorize(x.grains))))

ffq_hei_whi_c1 = pd.read_csv("../data/raw/whi/diet/ffq_hei_c1.txt", sep="\t",
                             skiprows=10)
ffq_hei_whi_c2 = pd.read_csv("../data/raw/whi/diet/ffq_hei_c2.txt", sep="\t",
                             skiprows=10)
ffq_hei_whi = (pd.concat([ffq_hei_whi_c1, ffq_hei_whi_c2], ignore_index=True)
               .rename({'SUBJID': 'subjID', 'VISITYR': 'visit_year', 'HEI2005': 'hei'}, 
                       axis=1)
               .filter(['subjID', 'visit_year', 'hei']))

diet_data_whi = ffq_nutrients_whi  # Room to merge in specific food intakes
diet_data_whi = (pd.merge(diet_data_whi, ffq_items_whi, on=["subjID",
                                                               "visit_year"]))
diet_data_whi = pd.merge(diet_data_whi, ffq_hei_whi, on=["subjID",
                                                         "visit_year"])

myhei_components = ["RM", "WG", "nuts", "FV", "fish", "dairy", "grains"]
diet_data_whi.loc[:, myhei_components] = (
    diet_data_whi.loc[:, myhei_components].apply(lambda x: x /
                                               diet_data_whi.tot_cal))
diet_data_whi = diet_data_whi.assign(myhei=lambda x: 3 * x.FV + x.fish + 2 * x.WG - x.grains + x.dairy + x.nuts - x.RM)

quant_diet_cols = ffq_nutrients_whi.columns
quant_diet_cols = diet_data_whi.columns[~diet_data_whi.columns.isin(["subjID",
                                                                     "visit_year"])]
diet_data_whi.loc[:, quant_diet_cols] = (
    diet_data_whi.loc[:, quant_diet_cols].apply(winsorize))

# Merge and save

whi_metadata = pd.merge(pheno_data_whi, diet_data_whi, 
                        on=["subjID", "visit_year"], how="inner")
whi_metadata.to_csv("../data/processed/metadata_whi.csv", index=False)

#GWAS prep

# Sample-to-subject mapping
sample_ids = (pd.read_csv("../data/raw/whi/gen/imputedSampleInfo.txt", sep="\t", skiprows=15)
              .rename({'SubjectID': 'subjID'}, axis=1)
              .filter(['SampleID', 'subjID']))

# Contains ancestry principal components linked to samples
sample_info_c1 = (pd.read_csv("../data/raw/whi/gen/sampleInfoWithPCs_c1.txt", sep="\t", skiprows=10)
                  .filter(['SAMPLE_ID','PC1','PC2','PC3','PC4','PC5']))
sample_info_c2 = (pd.read_csv("../data/raw/whi/gen/sampleInfoWithPCs_c2.txt", sep="\t", skiprows=10)
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


whi_metadata_medsAdj = (whi_metadata
                        .assign(ldl = lambda x: np.where(x.lipid_med, x.ldl /
                                                         0.75, x.ldl),
                                glu = lambda x: np.where(x.dm_med, x.glu /
                                                         0.75, x.glu),
                                sbp = lambda x: np.where(x.sbp, x.sbp  + 15,
                                                         x.sbp))
                        .copy())

#diet_vars = ['sfa', 'pufa', 'fat', 'f2c', 'palm']
#risk_factors = ["ldl", "tg", "glu", "sbp"]
#diet_vars = ["sfa", "f2c"]
#risk_factors = ["sbp", "bmi"]

diet_vars = ['sfa', 'pufa', 'fat', 'logf2c', 'palm', 'n3',
             'logsfa2pufa', 'logn62n3', 'myhei']
risk_factors = ["chol", "ldl", "hdl", "tg", "glu", "sbp", "bmi"]
#risk_factors.extend(["delta_" + rf for rf in risk_factors])

deltas = (whi_metadata_medsAdj
          .query('subjID in @pc_df.subjID')
          .filter(['subjID', 'visit_year'] + diet_vars + risk_factors)
          #.dropna(subset=diet_vars + risk_factors)
          .groupby('subjID')
          .filter(lambda x: x.shape[0] > 1)
          .sort_values('visit_year')
          .groupby('subjID')
          #.agg(lambda x: x.dropna().iloc[[0, -1]].diff().iloc[1] if x.count() >
          #     1 else np.nan)
          .agg(lambda x: x.iloc[-1] - x.iloc[0])
          .reset_index()
          .drop('visit_year', axis=1)
          .rename({var: "delta_" + var for var in diet_vars + risk_factors},
                  axis=1))

gwas_covars = ['dm_trial', 'sex', 'race', 'age',
               'PC1', 'PC2', 'PC3', 'PC4', 'PC5']
#gwas_covars.extend(diet_vars + risk_factors)

#delta_dv_rf_combos = [["delta_" + dv, "delta_" + rf] for dv, rf in dv_rf_combos]
#delta_dv_rf_combo_names = [dv + "_" + rf for dv, rf in delta_dv_rf_combos]

gwas_phenos = (whi_metadata_medsAdj
               .query('dm_trial == False')
               .query('visit_year == 0')
               .filter(['subjID'] + gwas_covars + diet_vars + risk_factors)
               #.dropna()
               .drop_duplicates(subset=["subjID"])
               .merge(deltas, on="subjID", how="left")
               .merge(pc_df, on="subjID")
               .assign(FID = lambda x: x.SampleID, IID = lambda x: x.SampleID)
               .filter(['FID', 'IID', 'sex'] + diet_vars + risk_factors + 
                       list(deltas.columns) + 
                       [v for v in gwas_covars if v != "sex"]))

diet_vars = diet_vars + ["delta_" + dv for dv in diet_vars]
risk_factors = risk_factors + ["delta_" + rf for rf in risk_factors]
dv_rf_combos = list(itertools.product(diet_vars, risk_factors))
dv_rf_combo_names = [dv + "_" + rf for dv, rf in dv_rf_combos]
for dv, rf in dv_rf_combos:  # Calculate scaled products for each diet/RF combo
    #gwas_phenos[dv + "_" + rf] = stats.mstats.winsorize(
    #    scale(gwas_phenos[dv]) * scale(gwas_phenos[rf]), limits=[0.01, 0.01])
    gwas_phenos[dv + "_" + rf] = scale(gwas_phenos[dv]) * scale(gwas_phenos[rf])

final_phenos = dv_rf_combo_names + diet_vars + risk_factors
gwas_phenos_INT = (gwas_phenos  # Inverse-normal transforms of phenotypes
                   .filter(final_phenos)
                   .apply(rank_INT)
                   .rename({p: p + "_INT" for p in final_phenos}, axis=1))
gwas_phenos = pd.concat([gwas_phenos, gwas_phenos_INT], axis=1) 


#gwas_covars = ['subjID', 'sex', 'age', 'race', 
#               'dm_trial', 'visit_year',
#               'PC1', 'PC2', 'PC3', 'PC4', 'PC5']
#input_phenos = ['sfa_pct', 'pufa_pct', 'fat_pct', 'cho', 'f2c', 'sbp',
#                'ldl']
#phenotypes = ['dsfa_dsbp', 'dsfa_dldl']
#
#deltas = (whi_metadata
#          .filter(['subjID', 'visit_year'] + input_phenos)
#          .dropna(subset=['sfa_pct', 'sbp'])
#          .groupby('subjID')
#          .filter(lambda x: x.shape[0] > 1)
#          .sort_values('visit_year')
#          .groupby('subjID')
#          .agg(lambda x: x.iloc[-1] - x.iloc[0])
#          .reset_index()
#          .drop('visit_year', axis=1)
#          .rename({'sfa_pct': 'delta_sfa_pct', 'sbp': 'delta_sbp', 'ldl':
#                   'delta_ldl'}, axis=1))
#
#gwas_phenos = (whi_metadata
#               .query('(ht_med == False | ht_med.isna()) & dm_trial == False')
#               .filter(gwas_covars)
#               .drop_duplicates(subset=['subjID'])
#               .merge(deltas, on="subjID")
#               .merge(pc_df, on="subjID")
#               .assign(dsfa_dsbp = lambda x: scale(x.delta_sfa_pct) *
#                       scale(x.delta_sbp),
#                       dsfa_dldl = lambda x: scale(x.delta_sfa_pct) *
#                       scale(x.delta_ldl))
#               .assign(FID = lambda x: x.SampleID, IID = lambda x: x.SampleID)
#               .filter(['FID', 'IID', 'sex'] + phenotypes +
#                       [v for v in gwas_covars if v != "sex"]))
#gwas_phenos_WIN = (gwas_phenos
#                   .filter(phenotypes)
#                   .apply(stats.mstats.winsorize, limits=[0.01, 0.01])
#                   .rename({p: p + "_WIN" for p in phenotypes}, axis=1))
#gwas_phenos_INT = (gwas_phenos
#                   .filter(phenotypes)
#                   .apply(rank_INT)
#                   .rename({p: p + "_INT" for p in phenotypes}, axis=1))
#gwas_phenos = pd.concat([gwas_phenos, gwas_phenos_WIN, gwas_phenos_INT], axis=1) 

gwas_phenos.query('race == "white"').to_csv(
    "../data/processed/whi/whi_white_gwas_phenos_long.txt", sep=" ", na_rep="NA", 
    index=False)
gwas_phenos.query('race == "black"').to_csv(
    "../data/processed/whi/whi_black_gwas_phenos_long.txt", sep=" ", na_rep="NA", 
    index=False)
gwas_phenos.query('race == "hispanic"').to_csv(
    "../data/processed/whi/whi_hispanic_gwas_phenos_long.txt", sep=" ", na_rep="NA", 
    index=False)

## List of WHI dietary modification trial partifcipants
#dm_phenos = (whi_metadata
#             .query('dm_trial == True')
#             .merge(pc_df, on="subjID")
#             .assign(FID = lambda x: x.SampleID,
#                     IID = lambda x: x.SampleID)
#             .filter(['FID', 'IID'])
#             .drop_duplicates())
#dm_phenos.to_csv("../data/processed/whi_DM_ids.txt", sep=" ", index=False)
