rf=$1
prefix=fat_pct_binary_${rf}_whi_white
all_prefix=${prefix}_all
nominal_prefix=${prefix}_nominalME
suggestive_prefix=${prefix}_suggestiveME
nominal_noSA_prefix=${prefix}_noSA_nominalME

#NOMINAL_SS=${NOMINAL_PREFIX}.res_annot

module load plink2

#plink2 --pfile ../whi/whi_DM \
#	--score $all_prefix.weights 1 6 9 \
#	--out dosage_scores/$all_prefix
#
#plink2 --pfile ../whi/whi_DM \
#	--score $nominal_prefix.weights 1 6 9 \
#	--out dosage_scores/$nominal_prefix
#
#plink2 --pfile ../whi/whi_DM \
#	--score $suggestive_prefix.weights 1 6 9 \
#	--out dosage_scores/$suggestive_prefix

ldpred_prefix=ldpred/ldpred_weights_ldl_LDpred

for f in -inf _p1.0000e+00 _p1.0000e-01 _p1.0000e-02 _p1.0000e-03 _p3.0000e-01 _p3.0000e-02 _p3.0000e-03; do
	plink2 --pfile ../whi/whi_DM \
		--score ${ldpred_prefix}${f}.txt 3 4 7 \
		--out dosage_scores/ldpred${f}_${rf}
done
