import sys
import pandas as pd
import statsmodels.formula.api as smf


rf, scorefile = sys.argv[1:]

dm_phenos = pd.read_csv("../data/processed/gen3/whi_white_DM_phenos.txt", sep=" ")

whi_sample_to_subject = (pd.read_csv("../data/raw/whi/sample_info.txt", sep="\t", skiprows=15)
			 .rename({'SampleID': 'sampleID', 'SubjectID': 'subjID'}, axis=1)
			 .filter(['sampleID', 'subjID']))

scores = (pd.read_csv(scorefile, delim_whitespace=True)
	  .rename({'IID': 'sampleID', 'SCORE': 'score'}, axis=1)
	  .merge(whi_sample_to_subject, on="sampleID")
	  .merge(dm_phenos, on="subjID"))

scores = scores.assign(delta_bmi = lambda x: x.delta_bmi.clip(lower=-5, upper=5),
		       delta_sbp = lambda x: x.delta_sbp.clip(lower=-30, upper=30),
		       delta_glu = lambda x: x.delta_glu.clip(lower=-40, upper=40),
		       delta_ldl = lambda x: x.delta_ldl.clip(lower=-60, upper=40),
		       delta_tg = lambda x: x.delta_tg.clip(lower=-75, upper=75))
print(smf.ols(f'delta_{rf} ~ score', data=scores.query('dm_intervention == True')).fit().summary().tables[1])
print(smf.ols(f'delta_{rf} ~ score', data=scores.query('dm_intervention == True & delta_sfa < delta_sfa.quantile(0.5)')).fit().summary().tables[1])
print(smf.ols(f'delta_{rf} ~ score + baseline_{rf}', data=scores.query('dm_intervention == True')).fit().summary().tables[1])
print(smf.ols(f'delta_{rf} ~ score * dm_intervention', data=scores).fit().summary().tables[1])
