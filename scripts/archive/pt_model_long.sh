#!/bin/bash


module load plink
module load plink2

DIR=$1
SUMSTATS=$DIR/all_meta.res
SNPANNO=../data/processed/snp_annotations/snp_annot_hg19_nodups.txt

## Generate bed/bim/fam to allow use of --clump
#plink2 --pfile ../data/processed/whi/whi \
#	--keep ../data/processed/whi_DM_ids.txt \
#	--hard-call-threshold 0.2 \
#	--make-bed \
#	--out ../data/processed/whi/whi_DM_tmp
#
#plink --bfile ../data/processed/whi/whi_DM_tmp \
#	--update-chr $SNPANNO 1 3 \
#	--update-map $SNPANNO 2 3 \
#	--make-bed \
#	--out ../data/processed/whi/whi_DM
#
#rm ../data/processed/whi/whi_DM_tmp*

# Same as above, but for WHI-only M-A results...
SUMSTATS_WHI=$DIR/whi_meta.res
plink --bfile ../data/processed/whi/whi_DM \
	--clump $SUMSTATS_WHI \
	--clump-field 'P' \
	--out $DIR/whi_score_weights

python << EOF
import pandas as pd
clumps = pd.read_csv("$DIR/whi_score_weights.clumped", delim_whitespace=True)
sum_stats = pd.read_csv("$SUMSTATS_WHI", delim_whitespace=True,
			usecols=['SNP', 'BETA'])
snp_annot = (pd.read_csv("../data/processed/snp_annotations/snp_annot_hg19_nodups.txt",
			sep="\t", header=None, names=['CHR', 'BP', 'SNP', 'REF', 'ALT', 'altname'])
	     .filter(['SNP', 'ALT']))
weights = (sum_stats[sum_stats.SNP.isin(clumps.SNP)]
 	   .merge(snp_annot, on="SNP")
	   .filter(['SNP', 'ALT', 'BETA']))
weights.to_csv("$DIR/whi_score_weights.txt", sep="\t", index=False)
EOF

# Calculate scores for WHI
plink2 --pfile ../data/processed/whi/whi \
    --score $DIR/whi_score_weights.txt 1 2 3 \
    --out $DIR/whi_scores_whi
