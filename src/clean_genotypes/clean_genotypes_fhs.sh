#!/bin/bash

module load gcta
module load plink
module load plink2
module load R

GENODIR=../data/processed/fhs

# Use dose2plink to convert minimac output to plink 1 dosage format
fhs_c1=/cluster/home/kweste01/kw/ncbi/dbGaP-60162_FHS_gen_IRBMDS/PhenoGenotypeFiles/ChildStudyConsentSet_phs000342.Framingham.v17.p10.c1.HMB-IRB-MDS/\
	GenotypeFiles/phg000679.v1.FHS_SHARe_Imputed_1000G.genotype-imputed-data.c1
fhs_c2=/cluster/home/kweste01/kw/ncbi/dbGaP-60163_FHS_gen_IRBNPUMDS/PhenoGenotypeFiles/ChildStudyConsentSet_phs000342.Framingham.v17.p10.c2.HMB-IRB-NPU-MDS/\
	GenotypeFiles/phg000679.v1.FHS_SHARe_Imputed_1000G.genotype-imputed-data.c2

#for ((i=1; i<=22; i++)); do
#	gunzip < $fhs_c1/imputed-metrics/machout.chr$i.info.gz > tmp_fhs_chr$i.info
#	awk '{if ($7 < 0.3) print $1}' tmp_fhs_chr$i.info > $GENODIR/fhs_chr${i}_low_qual_snps.txt
#	gunzip < $fhs_c1/machout.chr$i.dose_GRU.gz > tmp_fhs_chr${i}_c1.dose
#	gunzip < $fhs_c2/machout.chr$i.dose_NPU.gz > tmp_fhs_chr${i}_c2.dose
#	cat tmp_fhs_chr${i}_c1.dose tmp_fhs_chr${i}_c2.dose > tmp_fhs_chr$i.dose
#	awk '{print $1"->"$0}' tmp_fhs_chr$i.dose > tmp_fhs_chr$i.dosefix
#	~/kw/opt/dose2plink -dose tmp_fhs_chr$i.dosefix \
#		-info tmp_fhs_chr$i.info \
#		-out $GENODIR/fhs_chr$i
#	rm tmp_fhs*
#done

# "Custom" merge of chromosomes
#gunzip < $GENODIR/fhs_chr1.pdat.gz > $GENODIR/fhs.pdat
#for ((i=2; i<=22; i++)); do
#	gunzip < $GENODIR/fhs_chr$i.pdat.gz | tail -n +2 >> $GENODIR/fhs.pdat
#done

## Housekeeping prep for plink2 filters
#cat $GENODIR/fhs_chr*_low_qual_snps.txt > $GENODIR/fhs_low_qual_snps.txt
#awk -F, '{print $2,$2}' ../data/processed/metaData.csv > ../data/processed/meth_ids.txt
SNPANNO=../data/processed/snp_annotations/snp_annot_hg19_nodups.txt  # Contains 1000 Genomes SNP locations, rsIDs, and alleles from Ensembl

# Read-in, SNP annotations and R2 filtering, and pfile output using plink2
plink2 --import-dosage $GENODIR/fhs.pdat \
	--psam $GENODIR/fhs_chr1.pfam \
	--exclude $GENODIR/fhs_low_qual_snps.txt \
	--maf 0.001 \
	--update-name $SNPANNO 3 6 \
	--ref-allele $SNPANNO 4 3 \
	--make-pgen \
	--out $GENODIR/fhs

## SNP filtering by MAF and call rate and selection of IDs with methylation
#plink2 --pfile $GENODIR/fhs_unfiltered \
#	--maf 0.001 \
#	--geno 0.1 \
#	--make-pgen \
#	--out $GENODIR/fhs


#plink2 --pfile $GENODIR/fhs \
#	--make-bed \
#	--out $GENODIR/fhs

#tail -n +12 ../data/fhs/fhs_pedigree.txt | 
#	awk -F'\t' 'BEGIN{OFS=","} {print $4,$2,$5,$6,$9,-9}' |
#	sed 's/,,,/,0,0,/g' |  # Both parents missing
#	sed 's/,,/,0,/g' |  # One parent missing
#	tr ',' ' ' \
#	> ../data/processed/fhs/fhs_pedigree_plink.txt
#
#awk '{print $2,$2,$1,$2}' ../data/processed/fhs/fhs_pedigree_plink.txt > ../data/processed/fhs/fhs_id_update.txt
#plink --bfile $GENODIR/fhs \
#	--update-ids $GENODIR/fhs_id_update.txt \
#	--update-chr $SNPANNO 1 3 \
#	--update-map $SNPANNO 2 3 \
#	--make-bed \
#	--out $GENODIR/fhs

#awk '{print $1,$2,$3,$4}' ../data/processed/fhs/fhs_pedigree_plink.txt > ../int/fhs/fhs_parent_update.txt
#plink --bfile $GENODIR/fhs \
#	--update-parents $GENODIR/fhs_parent_update.txt \
#	--make-just-fam \
#	--out $GENODIR/fhs

#awk '!seen[$1]++' $GENODIR/fhs.fam > $GENODIR/unrelated_ids.txt

#plink --bfile $GENODIR/fhs \
#	#--filter-founders \
#	--keep $GENODIR/unrelated_ids.txt \
#	--make-bed \
#	--out $GENODIR/fhs_founders












# Loop through bfiles, removing SNPs with identical alleles and trimming to only SNPs w/ ACGT notation
#ls ../int/plink_imputed/*.bed | sed 's/.bed//g' > ../int/mergeList.txt
#mkdir ../int/plink_imputed_noIdentical
#for pref in $(cat ../int/mergeList.txt); do
#	prefOnly=$(echo $pref | cut -d '/' -f 4)
#	awk '$5==$6 {print $2}' ${pref}.bim > ../int/identicalAlleles.txt
#	plink --bfile $pref --exclude ../int/identicalAlleles.txt --alleleACGT --make-bed --out ../int/plink_imputed_noIdentical/tmp
#	plink --bfile ../int/plink_imputed_noIdentical/tmp --snps-only just-acgt --make-bed --out ../int/plink_imputed_noIdentical/${prefOnly}
#	rm ../int/plink_imputed_noIdentical/tmp.*
#done

# ls ../int/plink_imputed_noIdentical/*.bed | sed 's/.bed//g' > ../int/mergeList_noIdentical.txt
# #join -1 2 -2 1 -o 1.1,1.1, <(cat ${whi_c1_root}/phg000592.v1.WHI_Imputation.sample-info.MULTI/phg000592.v1_release_manifest.txt | tail -n +17 | cut -f 1,2 | sort -k 2,2) <(cat ../int/meth_ids.txt | sort -k 1,1) | awk '{print $1,$1}' > ../int/whi_sampIDs.txt
# plink --merge-list ../int/mergeList_noIdentical.txt --make-bed --out ../int/plink_imputed_noIdentical/noIdenticalMerge  # Merges all datasets into single plink set with only relevant subjects and performs variant QC
# mkdir ../int/plink_imputed_noMultiAllele
# for pref in $(cat ../int/mergeList_noIdentical.txt); do
#         newPref=../int/plink_imputed_noMultiAllele/$(echo $pref | cut -d '/' -f 4)
#         plink --bfile $pref --exclude ../int/plink_imputed_noIdentical/noIdenticalMerge-merge.missnp --make-bed --out $newPref
# 	cat ${newPref}.bim > tmpBim
# 	cut -f 2 ${newPref}.bim | awk -F ':' '{print $1,$2}' > tmpLoc
# 	paste <(awk '{print $1}' tmpLoc) <(cut -f 2 tmpBim) <(awk '{print 0,$2}' OFS='\t' tmpLoc) <(awk '{print $5,$6}' OFS='\t' tmpBim) > newBim && mv newBim ${newPref}.bim
# done
#ls ../int/plink_imputed_noMultiAllele/*.bed | sed 's/.bed//g' > ../int/mergeList_noMultiAllele.txt
#mkdir -p ../int/finalGenos
#cat ${whi_c1_root}/phg000592.v1.WHI_Imputation.sample-info.MULTI/phg000592.v1_release_manifest.txt | tail -n +17 | cut -f 1,2 | awk '{print $1,$1,$2,$2}' | sort -u -k 3 > ../int/idUpdate.txt
#plink --merge-list ../int/mergeList_noMultiAllele.txt --geno 0.01 --maf 0.01 --update-ids ../int/idUpdate.txt --indep 50 5 1.5 --make-bed --out ../int/finalGenos/unpruned  # Generates filtered plink set and list of SNPs to prune
#plink --bfile ../int/finalGenos/unpruned --extract ../int/finalGenos/unpruned.prune.in --recode A --out ../int/finalGenos/pruned  # Generates a .raw file with 0/1/2 genotype calls
