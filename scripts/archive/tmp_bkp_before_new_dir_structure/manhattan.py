import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from assocplots.manhattan import manhattan


ma_res_filename = sys.argv[1]

ma_res = (pd.read_csv(ma_res_filename, delim_whitespace=True)
          .query('P < 1e-2'))

cmap = plt.get_cmap('viridis')
colors = [cmap(i) for i in [0.0,0.33,0.67,0.90]]
chr_labels = np.array([str(i) for i in range(1, 23)])
chr_labels[12::2] = ''

manhattan(ma_res.P.values,
          ma_res.BP.values,
          ma_res.CHR.values,
          '',
          chrs_names=chr_labels,
          title="Meta-analysis results",
          colors=colors,
          lines=[5,7.5],
          lines_colors=['b','r'])
plt.savefig(ma_res_filename.replace(".meta", "_manhattan.png"))
          
