#!/bin/bash

module load plink
module load plink2

OUTDIR=$1

# LD-based clumping
plink --bfile ../int/plink_imputed/updated/whi_genos \  # LD estimates from WHI
       	--clump $OUTDIR/fhs_whi.meta \ 
	--clump-p1 0.0001 \
	--clump-r2 0.5 \
	--out $OUTDIR/fhs_whi

# Calculate clump-based scores
python << EOF
import pandas as pd
clumps = pd.read_csv("$OUTDIR/fhs_whi.clumped", delim_whitespace=True)
meta_res = pd.read_csv("$OUTDIR/fhs_whi.meta", delim_whitespace=True,
		       usecols=['SNP','A1','BETA'])
(meta_res[meta_res.SNP.isin(clumps.SNP)]
 .to_csv("$OUTDIR/clump_score_weights.txt", sep="\t", index=False)
EOF

calc_score () {
	local pop=$2
	local dir=$1
	plink --bfile bfile \
		--score $dir/clump_score_weights.txt 1 2 3 \
		--out $OUTDIR/$pop
}


calc_score ../int/plink_imputed/fhs_genos fhs
calc_score ../int/plink_imputed/updated/whi_genos whi
# also for bprhs once using dosages
# need score nuance w/ normalization once using dosages?

# Prep for model-based score development
echo -e "../int/plink_imputed/updated/whi_genos\n../int/plink_imputed/fhs_genos" > $OUTDIR/train_merge_list.txt
cat $OUTDIR/train_ids_whi.txt $OUTDIR/train_ids_fhs.txt > $OUTDIR/all_train_ids.txt
awk '{if ($7 < 0.001) print $3}' $OUTDIR/fhs_whi.meta > $OUTDIR/lasso_snps.txt
plink --merge-list $OUTDIR/train_merge_list.txt \
	--extract $OUTDIR/lasso_snps.txt \
	--keep $OUTDIR/all_train_ids.txt \
	--a2-allele ../int/snp_annotations/snp_annot_hg19_nodups.txt 4 3 \
	--make-bed \
	--recode AD vcf-iid \
	--out $OUTDIR/train_genos

# LASSO
plink --bfile train_genos \
	--lasso 0.2 \
	--out lasso_res

# Other model-based
# create separate python script to perform arbitrary machine learning modeling (training and score calculation?)
