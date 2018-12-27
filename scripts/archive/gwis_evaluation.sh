#!/bin/bash

module load python/3.6.0
module load plink/1.90b5

# Preliminaries
DIETVAR=$1
PHENO=$2
OUTDIR=$3
mkdir -p $OUTDIR
FHS_EXAM=$4

# Create train/test split in WHI
#awk -v trainfile=$OUTDIR/whi_nonDM_ids.txt -v testfile=$OUTDIR/whi_DM_ids.txt \
#	'{if (rand() < 0.7) {print $1,$2 > trainfile} else {print $1,$2 > testfile}}' \
#	../int/plink_imputed/updated/whi_genos.fam

#python << EOF
#import pandas as pd
#md_whi = (pd.read_csv("../int/metaData_whi.csv")
#	  .query('visitYear == 0'))
#md_whi.loc[lambda x: x.dm_trial == True, ['subjID', 'subjID']].to_csv(
#	"../int/whi_DM_ids.txt", sep=" ", header=False)
#md_whi.loc[lambda x: x.dm_trial == False, ['subjID', 'subjID']].to_csv(
#	"../int/whi_nonDM_ids.txt", sep=" ", header=False)
#EOF

# Prepare phenotypes
source ../env/bin/activate  # Turn on python venv
python clean_metadata_fhs.py $OUTDIR
python clean_metadata_whi.py $OUTDIR

# If run_gwis argument is given, submit jobarrays to run GWIS analysis in FHS & WHI
if [ "$5" == "run_gwis" ]
then
	# Run GWIS for each cohort/chromosome and concatenate chromosomes
	echo -n > $OUTDIR/fhs_complete.txt
	echo -n > $OUTDIR/whi_complete.txt
	sbatch --mem 10G -t 2:00:00 --array=1-22 --output=fhs%a.out run_gwis_chr_fhs.sh $DIETVAR $PHENO $OUTDIR $FHS_EXAM
	sbatch --mem 10G -t 2:00:00 --array=1-22 --output=whi%a.out run_gwis_chr_whi.sh $DIETVAR $PHENO $OUTDIR ancestry
	until [ `wc -l <$OUTDIR/fhs_complete.txt` == 22 ] && [ `wc -l <$OUTDIR/whi_complete.txt` == 22 ]; do sleep 10; done
	cat $OUTDIR/fhs_res_chr*_int.txt | awk 'NR==1 || !/^CHR/' > $OUTDIR/fhs_res_int.txt
	cat $OUTDIR/whi_res_chr*_int.txt | awk 'NR==1 || !/^CHR/' > $OUTDIR/whi_res_int.txt

	# Meta-analysis of cohort results
	plink --meta-analysis $OUTDIR/fhs_res_int.txt $OUTDIR/whi_res_int.txt + qt \
	--out $OUTDIR/fhs_whi

	# Clumping of meta-analysis hits (based on LD in full WHI)
	plink --bfile ../int/plink_imputed/updated/whi_genos \
	--clump $OUTDIR/fhs_whi.meta \
	--out $OUTDIR/fhs_whi
fi


# Generate basic stats on agreement of FHS and WHI GWIS results
python gwis_agreement.py $OUTDIR/fhs_res_int.txt $OUTDIR/whi_res_int.txt

# Generate a manhattan plot from M-A results
python manhattan.py $OUTDIR/fhs_whi.meta

# Python wrangling to generate an "interaction score" weights matrix
python << EOF
import pandas as pd
clumps = pd.read_csv("$OUTDIR/fhs_whi.clumped", delim_whitespace=True)
meta_res = pd.read_csv("$OUTDIR/fhs_whi.meta", delim_whitespace=True,
		       usecols=['SNP','A1','BETA'])
meta_res[meta_res.SNP.isin(clumps.SNP)].to_csv("$OUTDIR/score_weights.txt", sep="\t", index=False)
EOF

# Calculate interaction scores using plink
plink --bfile ../int/plink_imputed/updated/whi_genos \
	--score $OUTDIR/score_weights.txt \
	--out $OUTDIR/interaction_scores

# Run regressions to evaluate performance of interaction score
python test_interaction.py $DIETVAR $PHENO $OUTDIR

cd $OUTDIR
tar -cf ../../../output/${DIETVAR}_${PHENO}_res.tar agreement.txt fhs_whi_manhattan.png test_res.txt

deactivate  # Turn off python venv
