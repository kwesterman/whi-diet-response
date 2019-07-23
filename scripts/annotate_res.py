import sys
import pandas as pd


res_file = sys.argv[1]

res = (pd.read_csv(res_file, delim_whitespace=True)
	.dropna())

bim = pd.read_csv("../whi/whi_hardcalls.bim", delim_whitespace=True, header=None,
	names=['CHR', 'SNP', 'CM', 'BP', 'A1', 'A2'])

(pd.merge(res, bim, on="SNP", how="left", suffixes=['', '.bim'])
	.drop(['CHR.bim', 'CM', 'BP.bim', 'A1.bim'], axis=1)
	.to_csv(res_file + "_annot", sep="\t", index=False))
