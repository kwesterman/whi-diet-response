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
