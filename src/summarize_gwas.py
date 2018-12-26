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


res_file = sys.argv[1]


# Import results file
res = pd.read_csv(res_file, delim_whitespace=True)

# Assess inflation
def make_qqplot(pvec):
    lam = np.median(stats.chi2.ppf(1 - res.P, df=1)) / stats.chi2.ppf(0.5, df=1)
    x = -np.log10((np.arange(len(pvec)) + 1) / len(pvec))
    y = -np.log10(np.sort(pvec))
    plt.plot(x, x, color="r")
    plt.scatter(x, y, color="blue", s=0.1)
    plt.xlabel("-logP expected")
    plt.ylabel("-logP observed")
    plt.title(f"lambda = {lam}")
    return None

make_qqplot(res.P)
plt.savefig(res_file.replace(".res", "_qq.png"))

#qqplot([res.P],
#       ['GWAS results'],
#       color=['b'],
#       title='')
#plt.savefig(res_filename.replace(".txt", "_qqplot.png"))

#sm.qqplot(res.P, stats.uniform, line='45', ms=0.5)
#plt.savefig(res_file.replace(".res", "_qq.png"))

# SNP annotations
snp_annot = pd.read_csv(
    "../data/processed/snp_annotations/snp_annot_hg19_nodups.txt",
    sep="\t", header=None, usecols=[0, 1, 2], names=['chr', 'bp', 'id'])

res = (res 
       .query('P < 1e-2') 
       .merge(snp_annot, left_on="ID", right_on="id"))

# Set color scheme and create Manhattan plot
cmap = plt.get_cmap('viridis')
colors = [cmap(i) for i in [0.0,0.33,0.67,0.90]]
chr_labels = np.array([str(i) for i in range(1, 23)])
chr_labels[12::2] = ''

manhattan(res.P.values,
          res.bp.values,
          res.chr.astype(str),
          '',
          chrs_names=chr_labels,
          title='',
          colors=colors,
          lines=[5,7.5],
          lines_colors=['b','r'])
plt.savefig(res_file.replace(".res", "_manhattan.png"))
