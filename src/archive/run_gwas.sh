#!/bin/bash

module load python/3.6.0
module load plink/1.90b5

# Preliminaries
POP=$1

DIR=../int/sfa

echo -n > $DIR/${POP}_complete.txt
sbatch --mem 10G -t 2:00:00 --array=1-22 --output=${POP}%a.out run_gwas/run_gwas_chr_$POP.sh $DIR
until [ `wc -l <$DIR/${POP}_complete.txt` == 22 ]; do sleep 10; done
cat $DIR/${POP}_res_chr*.txt | awk 'NR==1 || !/^CHR/' > $DIR/${POP}_res.txt
#plink --bfile $GENODIR/${POP}_genos \
#	--keep $DIR/${POP}_phenos.txt \
#	--clump $DIR/${POP}_res.txt \
#	--out $DIR/${POP}

## Python wrangling to generate an "interaction score" weights matrix
#python << EOF
#import pandas as pd
#clumps = pd.read_csv("$OUTDIR/fhs_whi.clumped", delim_whitespace=True)
#meta_res = pd.read_csv("$OUTDIR/fhs_whi.meta", delim_whitespace=True,
#		       usecols=['SNP','A1','BETA'])
#meta_res[meta_res.SNP.isin(clumps.SNP)].to_csv("$OUTDIR/score_weights.txt", sep="\t", index=False)
#EOF

## Calculate interaction scores using plink
#plink --bfile ../int/plink_imputed/updated/whi_genos \
#	--score $OUTDIR/score_weights.txt \
#	--out $OUTDIR/interaction_scores
