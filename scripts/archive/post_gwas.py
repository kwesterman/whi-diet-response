import sys
import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from scipy import stats
import statsmodels.api as sm
from assocplots.manhattan import manhattan
from assocplots.qqplot import *
import re


res_file = sys.argv[1]

# SNP annotations
snp_annot = (pd.read_pickle("../data/processed/snp_annotations/snp_annot_hg19_nodups.pkl")
             .set_index("id"))
#snp_annot = pd.read_csv(
#    "../data/processed/snp_annotations/snp_annot_hg19_nodups.txt",
#    sep="\t", header=None, usecols=[0, 1, 2], names=['chr', 'bp', 'id'])

# Import results file and merge with SNP annotation
res = (pd.read_csv(res_file, delim_whitespace=True)
       .rename({'ID': 'SNP'}, axis=1))  # GWAS and M-A have diff. ID col names
       #.merge(snp_annot, left_on="ID", right_on="id")
       #.loc[lambda x: x.chr.isin(np.arange(1, 23))])
res_annot = (res.set_index("SNP")
             .join(snp_annot, how="inner")
             .rename_axis("SNP")
             .reset_index()
             .sort_values(['chr', 'bp']))

# Meta-analysis-specific steps
if "meta" in res_file:
    full_N = np.max(res_annot.N)
    res_annot = (res_annot
                 .query('N == @full_N')
                 .filter(["SNP", "chr", "bp", "REF", "ALT", 
			  "N", "P", "P(R)", "BETA", "BETA(R)", "Q", "I"]))
    meta_res = res_annot.rename({'chr': 'CHR', 'bp': 'BP'}, axis=1)
    meta_res.to_csv(res_file.replace(".meta", ".res"), sep="\t", index=False)
    #res_annot['P'] = res_annot['P(R)']  # If want plots to reflect RE p-vals

# Assess inflation
def make_qqplot(pvec):
    lam = np.nanmedian(stats.chi2.ppf(1 - pvec, df=1)) / stats.chi2.ppf(0.5, df=1)
    x = -np.log10((np.arange(len(pvec)) + 1) / len(pvec))
    y = -np.log10(np.sort(pvec))
    plt.plot(x, x, color="r")
    plt.scatter(x, y, color="blue", s=0.1)
    plt.xlabel("-logP expected")
    plt.ylabel("-logP observed")
    plt.title(f"lambda = {lam.round(2)}")
    return None

make_qqplot(res_annot["P"])
plt.savefig(re.sub(r"\.[a-z]+", "_qq.png", res_file))

# Set color scheme and create Manhattan plot
res_annot = res_annot.loc[res_annot["P"] < 1e-2]  # Easier to plot fewer points

cmap = plt.get_cmap('viridis')
colors = [cmap(i) for i in [0.0,0.33,0.67,0.90]]
chr_labels = np.array([str(i) for i in range(1, 23)])
chr_labels[12::2] = ''

manhattan(res_annot["P"].values,
          res_annot.bp.values,
          res_annot.chr.astype(str),
          '',
          chrs_names=chr_labels,
          title='',
          colors=colors,
          lines=[5,7.5],
          lines_colors=['b','r'])
plt.savefig(re.sub(r"\.[a-z]+", "_manhattan.png", res_file))
