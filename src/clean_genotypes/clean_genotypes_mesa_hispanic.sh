#!/bin/bash

module load plink
module load plink2

DATADIR=/cluster/home/kweste01/kw/ncbi/mesa/MESA_SHARe_imputation_HIS_1kg_08232012
DESTDIR=/cluster/home/kweste01/kw/diet_response/data/processed/mesa

mkdir -p $DESTDIR

cat <(echo -e "ID_1 ID_2 missing\n0 0 0") \
	<(awk '{print $1,$1,0}' $DATADIR/imputation_sample_file_gen3_HIS.txt) \
	> $DESTDIR/mesa_hispanic.sample
cat $DATADIR/imputation_data/imp2_HIS_chr*_ref1kg > $DESTDIR/mesa_hispanic.gen

cat $DATADIR/SNP_info/imp2_HIS_chr*_ref1kg_info > $DESTDIR/mesa_hispanic.snpinfo
awk '$5>0.3 {print $2}' $DESTDIR/mesa_hispanic.snpinfo > $DESTDIR/mesa_hispanic_good_snps.txt

awk -F, '{print $1,$1}' ../data/processed/metadata_mesa.csv > $DESTDIR/mesa_ids.txt

#gunzip < $DESTDIR/mesa_hispanic.gen.gz > $DESTDIR/mesa_hispanic.gen
plink2 --gen $DESTDIR/mesa_hispanic.gen \
	--allow-extra-chr \
	--sample $DESTDIR/mesa_hispanic.sample \
	--extract $DESTDIR/mesa_hispanic_good_snps.txt \
	--keep $DESTDIR/mesa_ids.txt \
	--sort-vars \
	--make-pgen \
	--out $DESTDIR/mesa_hispanic_tmp
#rm $DESTDIR/mesa_hispanic.gen

#awk '{if($1=="X" || $1=="---") print $3}' $DESTDIR/mesa_hispanic.pvar > $DESTDIR/mesa_hispanic.rmsnps
mv $DESTDIR/mesa_hispanic_tmp.pvar $DESTDIR/mesa_hispanic_tmp.pvartmp
echo "#CHROM  POS     ID      REF     ALT" > $DESTDIR/mesa_hispanic_tmp.pvar
tail -n +2 $DESTDIR/mesa_hispanic_tmp.pvartmp | awk -v OFS='\t' '{print 0,0,$3,$4,$5}' \
	>> $DESTDIR/mesa_hispanic_tmp.pvar

SNPANNO=../data/processed/snp_annotations/snp_annot_hg19_nodups.txt
plink2 --pfile $DESTDIR/mesa_hispanic_tmp \
	--maf 0.001 \
	--make-pgen \
	--out $DESTDIR/mesa_hispanic
	#--allow-extra-chr \
	#--update-chr $SNPANNO 1 3 \
	#--update-map $SNPANNO 2 3 \

#join -1 2 -2 3 -o 2.1,2.3,2.2,1.4,1.5,1.6 \
#	<(sort -k2,2 $DESTDIR/mesa_hispanic.bim) \
#	<(sort -k3,3 ../data/processed/snp_annotations/snp_annot_hg19_nodups.txt) \
#	> $DESTDIR/mesa_hispanic.bim

#echo "$(awk '{print 24,$2,$3,$4,$5,$6}' $DESTDIR/mesa_hispanic.bim)" > $DESTDIR/mesa_hispanic.bim
## Replace "---" chromosome numbers with something plink likes better (to be overwritten below)
#
#plink --bfile $DESTDIR/mesa_hispanic \
#	--update-chr ../data/processed/snp_annotations/snp_annot_hg19_nodups.txt 1 3 \
#	--update-map ../data/processed/snp_annotations/snp_annot_hg19_nodups.txt 2 3 \
#	--maf 0.001 \
#	--make-bed \
#	--out $DESTDIR/mesa_hispanic
#
#plink --bfile $DESTDIR/mesa_hispanic \
#	--chr 1-22 \
#	--make-bed \
#	--out $DESTDIR/mesa_hispanic
