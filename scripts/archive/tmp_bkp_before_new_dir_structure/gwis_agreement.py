import sys
import pandas as pd
from scipy.stats import pearsonr


fhs_int_res_file = sys.argv[1]
whi_int_res_file = sys.argv[2]

print("Reading interaction result files...")
fhs_int_res = (pd.read_csv(fhs_int_res_file, 
                           delim_whitespace=True, 
                           usecols=['SNP','BETA','P'])
               .assign(BETA=lambda x: pd.to_numeric(x.BETA, errors="coerce"),
                       P=lambda x: pd.to_numeric(x.P, errors="coerce")))
whi_int_res = (pd.read_csv(whi_int_res_file, 
                           delim_whitespace=True, 
                           usecols=['SNP','BETA','P'])
               .assign(BETA=lambda x: pd.to_numeric(x.BETA, errors="coerce"),
                       P=lambda x: pd.to_numeric(x.P, errors="coerce")))

print("Merging...")
all_int_res = (pd.merge(fhs_int_res, whi_int_res,
                        on="SNP", suffixes=("_fhs","_whi"))
               .dropna())

trimmed = all_int_res.query('P_fhs<0.01 and P_whi<0.01')

num_shared_snps = all_int_res.shape[0]
num_shared_snps_p01 = trimmed.shape[0]

beta_corr = pearsonr(trimmed["BETA_fhs"], trimmed["BETA_whi"])

print("Writing results...")
with open(fhs_int_res_file.replace(
    "fhs_res_int.txt", "agreement.txt"), "w") as f:
    f.write(f'{num_shared_snps} SNPs overlapping.\n')
    f.write(f'{num_shared_snps_p01} SNPs overlapping at p < 0.01.\n')
    f.write(f'Correlation between beta coefficients with p<0.01 in both cohorts: {round(beta_corr[0],4)}.')
