#!/bin/bash

# This script assumes we start with available .pdat and .pfam files (plink1 dosage files) for each of the groups and merges them into a single plink2 pfile set while retaining dosages

module load plink
module load plink2

DIR=../data/processed/whi

declare -a groups=("as264" "garnet" "gecco_cyto" "gecco_init" "hipfx" "whims" "share_aa" "share_ha")

echo -n > $DIR/whi_combined.pfam
tail -n +2 $DIR/whi_as264.pdat | awk '{print $1,$2,$3}' | sort \
	> $DIR/shared_snps.txt

# Loop through groups to prep for later merge
for group in ${groups[@]}; do
	echo $group'...'
	awk '{print $1,$2,$3}' $DIR/whi_$group.pdat | sort \
		> $DIR/whi_$group.snpset  # collects SNP IDs and alleles
	awk '{print 0,$1,0,0}' $DIR/whi_$group.pdat \
		> $DIR/whi_$group.map
	cat $DIR/whi_$group.pfam >> $DIR/whi_combined.pfam
	echo "$(comm -12 $DIR/shared_snps.txt $DIR/whi_$group.snpset)" \
		> $DIR/shared_snps.txt
done

echo "$(cut -d ' ' -f 1 $DIR/shared_snps.txt)" > $DIR/shared_snps.txt  # SNPs that are shared AND have same alleles
echo -n > $DIR/dosage_files.txt

# Extract shared SNPs with matching alleles for all groups
for group in ${groups[@]}; do
	plink --dosage $DIR/whi_$group.pdat format=1 \
		--fam $DIR/whi_$group.pfam \
		--map $DIR/whi_$group.map \
		--extract $DIR/shared_snps.txt \
		--write-dosage \
		--out $DIR/whi_$group
	echo $DIR'/whi_'$group'.out.dosage' >> $DIR/dosage_files.txt
done

# Merge all dosage files
plink --dosage $DIR/dosage_files.txt list format=1 \
	--fam $DIR/whi_combined.pfam \
	--write-dosage \
	--out $DIR/whi

cat $DIR/whi_*_low_qual_snps.txt > $DIR/whi_low_qual_snps.txt

# Import dosage and general cleanup using plink2
SNPANNO=../data/processed/snp_annotations/snp_annot_hg19_nodups.txt
plink2 --import-dosage $DIR/whi.out.dosage \
	--psam $DIR/whi_combined.pfam \
	--exclude $DIR/whi_low_qual_snps.txt \
	--update-name $SNPANNO 3 6 \
	--ref-allele $SNPANNO 4 3 \
	--maf 0.001 \
	--make-pgen \
	--out $DIR/whi_tmp

# while --alleleACGT is not yet implemented in plink2...
mv $DIR/whi_tmp.pvar $DIR/whi_tmp.pvartmp
awk '{gsub(1,"A",$4); gsub(2,"C",$4); gsub(3,"G",$4); gsub(4,"T",$4); \
	gsub(1,"A",$5); gsub(2,"C",$5); gsub(3,"G",$5); gsub(4,"T",$5); print}' \
	< $DIR/whi_tmp.pvartmp \
	> $DIR/whi_tmp.pvar

plink2 --pfile $DIR/whi_tmp \
	--ref-allele $SNPANNO 4 3 \
	--make-pgen \
	--out $DIR/whi
