import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from assocplots.manhattan import manhattan


res_filename = sys.argv[1]
plot_title = sys.argv[2]

res = (pd.read_csv(res_filename, delim_whitespace=True)
       .assign(P=lambda x: pd.to_numeric(x.P, errors="coerce"),
               BP=lambda x: pd.to_numeric(x.BP, errors="coerce"))
       .query('P < 1e-2'))

cmap = plt.get_cmap('viridis')
colors = [cmap(i) for i in [0.0,0.33,0.67,0.90]]
chr_labels = np.array([str(i) for i in range(1, 23)])
chr_labels[12::2] = ''

manhattan(res.P.values,
          res.BP.values,
          res.CHR.astype(str),
          '',
          chrs_names=chr_labels,
          title=plot_title,
          colors=colors,
          lines=[5,7.5],
          lines_colors=['b','r'])
plt.savefig(res_filename.replace(".txt", "_manhattan.png"))
          
