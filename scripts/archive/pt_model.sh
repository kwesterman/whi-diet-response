#!/bin/bash

module load python/3.6.0
module load plink/1.90b5

# Preliminaries
SUMSTATS=$1
PLINKSET=$2
OUT_PREFIX=$3

# Run clumping procedure (pruning + thresholding)
plink --bfile $PLINKSET \
	--clump $SUMSTATS \
	--out $OUT_PREFIX

# Python wrangling to generate a P/T weights file
python << EOF
import pandas as pd
clumps = pd.read_csv("$OUT_PREFIX.clumped", delim_whitespace=True)
sum_stats = pd.read_csv("$SUMSTATS", delim_whitespace=True,
			usecols=['SNP', 'A1', 'BETA'])
weights = sum_stats[sum_stats.SNP.isin(clumps.SNP)]
weights.to_csv("${OUT_PREFIX}_weights.txt", sep="\t", index=False)
EOF
