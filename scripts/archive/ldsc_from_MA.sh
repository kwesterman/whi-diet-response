#!/bin/bash


module load anaconda

LDSC_DIR=~/kw/opt/ldsc
LDSC_REF_DIR=../data/processed/architecture/ldsc_gen_corr
SUMSTATS=$1

source activate ~/kw/opt/conda-envs/ldsc

python $LDSC_DIR/munge_sumstats.py \
	--sumstats $SUMSTATS \
	--N 5000 \
	--a1 ALT \
	--a2 REF \
	--merge-alleles $LDSC_REF_DIR/w_hm3.snplist \
	--out $SUMSTATS

python $LDSC_DIR/ldsc.py \
	--h2 $SUMSTATS.sumstats.gz \
	--ref-ld-chr $LDSC_REF_DIR/eur_w_ld_chr/ \
	--w-ld-chr $LDSC_REF_DIR/eur_w_ld_chr/ \
	--out $SUMSTATS.ldsc
