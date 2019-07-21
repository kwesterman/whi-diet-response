#!/bin/bash

# This script assumes we start with available .pdat and .pfam files (plink1 dosage files) for each of the groups and merges them into a single plink2 pfile set while retaining dosages

module load plink
module load plink2

DIR=../data/processed/whi

declare -a groups=("as264" "garnet" "gecco_cyto" "gecco_init" "hipfx" "whims" "share_aa" "share_ha")

#echo -n > $DIR/whi_combined.pfam
#tail -n +2 $DIR/whi_${groups[0]}.pdat | awk '{print $1,$2,$3}' | sort \
#	> $DIR/shared_snps.txt
echo -n > $DIR/whi_bfiles.txt

# Loop through groups to prep for later merge
for group in ${groups[@]}; do
	echo $group'...'
	#awk '{print $1,$2,$3}' $DIR/whi_$group.pdat | sort \
	#	> $DIR/whi_$group.snpset  # collects SNP IDs and alleles
	#awk '{print 0,$1,0,0}' $DIR/whi_$group.pdat \
	#	> $DIR/whi_$group.map
done

SNPANNO=../data/processed/snp_annotations/snp_annot_hg19_nodups.txt

# Extract shared SNPs with matching alleles for all groups
for group in ${groups[@]}; do
	#plink --dosage $DIR/whi_$group.pdat format=1 \
	#	--fam $DIR/whi_$group.pfam \
	#	--map $DIR/whi_$group.map \
	#	--exclude $DIR/whi_${group}_low_qual_snps.txt \
	#	--allow-no-sex \
	#	--write-dosage \
	#	--out $DIR/whi_$group
	#	#--hard-call-threshold 0.1 \
	#	#--update-chr $SNPANNO 1 6 \
	#	#--update-map $SNPANNO 2 6 \
	#	#--chr 1-22 \


	#plink2 --import-dosage $DIR/whi_${group}.out.dosage \
	#	--psam $DIR/whi_${group}.pfam \
	#	--update-name $SNPANNO 3 6 \
	#	--export vcf vcf-dosage=DS \
	#	--out $DIR/whi_${group}

	#plink2 --import-dosage $DIR/whi_${group}.out.dosage \
	#	--psam $DIR/whi_${group}.pfam \
	#	--update-name $SNPANNO 3 6 \
	#	--hard-call-threshold 0.1 \
	#	--make-bed \
	#	--out $DIR/whi_${group}

	#plink --bfile $DIR/whi_${group} \
	#	--update-chr $SNPANNO 1 3 \
	#	--update-map $SNPANNO 2 3 \
	#	--alleleACGT \
	#	--make-bed \
	#	--out $DIR/whi_${group}

	echo "$DIR/whi_$group" >> $DIR/whi_bfiles.txt
done

#plink --merge-list $DIR/whi_bfiles.txt \
#	--a1-allele $SNPANNO 4 3 \
#	--make-bed \
#	--out $DIR/whi_hardcalls

#echo -n > $DIR/tmp_whi_bfiles.txt
#for group in ${groups[@]}; do
#	plink --bfile $DIR/whi_$group \
#		--exclude $DIR/whi_hardcalls-merge.missnp \
#		--make-bed \
#		--out $DIR/tmp_whi_${group}
#
#	echo "$DIR/tmp_whi_$group" >> $DIR/tmp_whi_bfiles.txt
#done

plink --merge-list $DIR/tmp_whi_bfiles.txt \
	--a2-allele $SNPANNO 4 3 \
	--make-bed \
	--out $DIR/whi_hardcalls

rm $DIR/tmp_*
#plink --bfile whi_test_merge \
#	--update-chr $SNPANNO 1 3 \
#	--update-map $SNPANNO 2 3 \
#	--alleleACGT \
#	--a2-allele $SNPANNO 4 3 \
#	--make-bed \
#	--out $DIR/whi_test_merge_annot






#echo "$(cut -d ' ' -f 1 $DIR/shared_snps.txt)" > $DIR/shared_snps.txt  # SNPs that are shared AND have same alleles
#echo -n > $DIR/dosage_files.txt
#
## Merge all dosage files
#plink --dosage $DIR/dosage_files.txt list format=1 \
#	--fam $DIR/whi_combined.pfam \
#	--write-dosage \
#	--out $DIR/whi
#
#cat $DIR/whi_*_low_qual_snps.txt > $DIR/whi_low_qual_snps.txt
#
## Import dosage and general cleanup using plink2
#SNPANNO=../data/processed/snp_annotations/snp_annot_hg19_nodups.txt
#plink2 --import-dosage $DIR/whi.out.dosage \
#	--psam $DIR/whi_combined.pfam \
#	--exclude $DIR/whi_low_qual_snps.txt \
#	--update-name $SNPANNO 3 6 \
#	--ref-allele $SNPANNO 4 3 \
#	--maf 0.001 \
#	--make-pgen \
#	--out $DIR/whi_tmp
#
## while --alleleACGT is not yet implemented in plink2...
#mv $DIR/whi_tmp.pvar $DIR/whi_tmp.pvartmp
#awk '{gsub(1,"A",$4); gsub(2,"C",$4); gsub(3,"G",$4); gsub(4,"T",$4); \
#	gsub(1,"A",$5); gsub(2,"C",$5); gsub(3,"G",$5); gsub(4,"T",$5); print}' \
#	< $DIR/whi_tmp.pvartmp \
#	> $DIR/whi_tmp.pvar
#
#plink2 --pfile $DIR/whi_tmp \
#	--ref-allele $SNPANNO 4 3 \
#	--make-pgen \
#	--out $DIR/whi
#
## If desired, some further cleaning/filtering and conversion to bfiles for use by other software
#plink2 --pfile $DIR/whi \
#	--make-bed \
#	--hard-call-threshold 0.1 \
#	--out $DIR/whi
#
#plink --bfile $DIR/whi \
#	--update-chr $SNPANNO 1 3 \
#	--update-map $SNPANNO 2 3 \
#	--make-bed \
#	--out $DIR/whi
#
#plink --bfile $DIR/whi \
#	--maf 0.01 \
#	--geno 0.1 \
#	--chr 1-22 \
#	--make-bed \
#	--out $DIR/whi
