import sys
import pandas as pd


filename = sys.argv[1]

res_annot = pd.read_csv(filename, sep="\t")
res_annot.dropna().to_csv("ldpred/" + filename + ".ldpred_input", sep="\t", index=False)
