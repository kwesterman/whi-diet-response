#!/bin/bash


module load plink
module load plink2

DIR=$1
RES_FILE=$2
#SUMSTATS=$DIR/all_meta.res
SNPANNO=../data/processed/snp_annotations/snp_annot_hg19_nodups.txt

## Generate bed/bim/fam to allow use of --clump
#plink2 --pfile ../data/processed/whi/whi \
#	--keep ../data/processed/whi_white_DM_ids.txt \
#	--hard-call-threshold 0.2 \
#	--make-bed \
#	--out ../data/processed/whi/whi_white_DM_tmp
#
#plink --bfile ../data/processed/whi/whi_white_DM_tmp \
#	--update-chr $SNPANNO 1 3 \
#	--update-map $SNPANNO 2 3 \
#	--make-bed \
#	--out ../data/processed/whi/whi_white_DM
#
#rm ../data/processed/whi/whi_white_DM_tmp*


# White results
SUMSTATS_WHITE=$DIR/$RES_FILE
plink --bfile ../data/processed/whi/whi_white_DM \
	--clump $SUMSTATS_WHITE \
	--clump-snp-field 'SNP' 'ID' \
	--clump-field 'P' \
	--out $DIR/white_score_weights

python << EOF
import pandas as pd
clumps = pd.read_csv("$DIR/white_score_weights.clumped", delim_whitespace=True)
sum_stats = pd.read_csv("$SUMSTATS_WHITE", delim_whitespace=True,
			usecols=['SNP', 'BETA'])
snp_annot = (pd.read_csv("../data/processed/snp_annotations/snp_annot_hg19_nodups.txt",
			sep="\t", header=None, names=['CHR', 'BP', 'SNP', 'REF', 'ALT', 'altname'])
	     .filter(['SNP', 'ALT']))
weights = (sum_stats[sum_stats.SNP.isin(clumps.SNP)]
 	   .merge(snp_annot, on="SNP")
	   .filter(['SNP', 'ALT', 'BETA']))
weights.to_csv("$DIR/white_score_weights.txt", sep="\t", index=False)
EOF

# Calculate white M-A scores for WHI
plink2 --pfile ../data/processed/whi/whi \
    --score $DIR/white_score_weights.txt 1 2 3 \
    --out $DIR/whi_scores_white

#### NOW REPEAT BUT WITH RANDOM EFFECTS ESTIMATES
## Use clump to generate score weights
#plink --bfile ../data/processed/whi/whi_DM \
#	--clump $SUMSTATS \
#	--clump-field 'P(R)' \
#	--out $DIR/score_weights_RE
#
#
## Python wrangling to generate a P/T weights file
#python << EOF
#import pandas as pd
#clumps = pd.read_csv("$DIR/score_weights_RE.clumped", delim_whitespace=True)
#sum_stats = pd.read_csv("$SUMSTATS", delim_whitespace=True,
#			usecols=['SNP', 'BETA(R)'])
#snp_annot = (pd.read_csv("../data/processed/snp_annotations/snp_annot_hg19_nodups.txt",
#			sep="\t", header=None, names=['CHR', 'BP', 'SNP', 'REF', 'ALT', 'altname'])
#	     .filter(['SNP', 'ALT']))
#weights = (sum_stats[sum_stats.SNP.isin(clumps.SNP)]
# 	   .merge(snp_annot, on="SNP")
#	   .filter(['SNP', 'ALT', 'BETA(R)']))
#weights.to_csv("$DIR/score_weights_RE.txt", sep="\t", index=False)
#EOF
#
## Calculate scores for WHI
#plink2 --pfile ../data/processed/whi/whi \
#    --score $DIR/score_weights_RE.txt 1 2 3 \
#    --out $DIR/whi_scores_RE

