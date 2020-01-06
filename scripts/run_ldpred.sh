#!/bin/bash


RF=$1
GENOS=../whi/whi_white_DM
SSFILE=ldpred/fat_pct_binary_${RF}_whi_white.res_annot.ldpred_input

#ldpred coord \
#	--gf $GENOS \
#	--ssf $SSFILE \
#	--out ldpred/ldpred_coord_${RF} \
#	--only-hm3 \
#	--eff_type LINREG \
#	--maf 0.000001 \
#	--rs SNP \
#	--A1 ALT \
#	--A2 REF \
#	--pos bp \
#	--chr chr \
#	--pval P \
#	--eff BETA \
#	--ncol OBS_CT
#	#--vgf $GENOS \

ldpred gibbs \
	--cf ldpred/ldpred_coord_${RF} \
	--ldr 500 \
	--ldf ldpred/whi_white_DM_${RF} \
	--out ldpred/ldpred_weights_${RF} \
	--use-gw-h2
