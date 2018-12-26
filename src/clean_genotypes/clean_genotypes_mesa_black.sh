#!/bin/bash

module load plink
module load plink2

DATADIR=/cluster/home/kweste01/kw/ncbi/mesa/MESA_SHARe_imputation_AFA_1kg_11202012
DESTDIR=/cluster/home/kweste01/kw/diet_response/data/processed/mesa

mkdir -p $DESTDIR

cat <(echo -e "ID_1 ID_2 missing\n0 0 0") \
	<(awk '{print $1,$1,0}' $DATADIR/imputation_sample_file_gen3_AFA.txt) \
	> $DESTDIR/mesa_black.sample
cat $DATADIR/imputation_data/imp2_AFA_chr*_ref1kg > $DESTDIR/mesa_black.gen

cat $DATADIR/SNP_info/imp2_AFA_chr*_ref1kg_info > $DESTDIR/mesa_black.snpinfo
awk '$5>0.3 {print $2}' $DESTDIR/mesa_black.snpinfo > $DESTDIR/mesa_black_good_snps.txt

awk -F, '{print $1,$1}' ../data/processed/metadata_mesa.csv > $DESTDIR/mesa_ids.txt

#gunzip < $DESTDIR/mesa_black.gen.gz > $DESTDIR/mesa_black.gen
plink2 --gen $DESTDIR/mesa_black.gen \
	--allow-extra-chr \
	--sample $DESTDIR/mesa_black.sample \
	--extract $DESTDIR/mesa_black_good_snps.txt \
	--keep $DESTDIR/mesa_ids.txt \
	--sort-vars \
	--make-pgen \
	--out $DESTDIR/mesa_black_tmp
#rm $DESTDIR/mesa_black.gen

#awk '{if($1=="X" || $1=="---") print $3}' $DESTDIR/mesa_black.pvar > $DESTDIR/mesa_black.rmsnps
mv $DESTDIR/mesa_black_tmp.pvar $DESTDIR/mesa_black_tmp.pvartmp
echo "#CHROM  POS     ID      REF     ALT" > $DESTDIR/mesa_black_tmp.pvar
tail -n +2 $DESTDIR/mesa_black_tmp.pvartmp | awk -v OFS='\t' '{print 0,0,$3,$4,$5}' \
	>> $DESTDIR/mesa_black_tmp.pvar

SNPANNO=../data/processed/snp_annotations/snp_annot_hg19_nodups.txt
plink2 --pfile $DESTDIR/mesa_black_tmp \
	--maf 0.001 \
	--make-pgen \
	--out $DESTDIR/mesa_black
	#--allow-extra-chr \
	#--update-chr $SNPANNO 1 3 \
	#--update-map $SNPANNO 2 3 \

#join -1 2 -2 3 -o 2.1,2.3,2.2,1.4,1.5,1.6 \
#	<(sort -k2,2 $DESTDIR/mesa_black.bim) \
#	<(sort -k3,3 ../data/processed/snp_annotations/snp_annot_hg19_nodups.txt) \
#	> $DESTDIR/mesa_black.bim

#echo "$(awk '{print 24,$2,$3,$4,$5,$6}' $DESTDIR/mesa_black.bim)" > $DESTDIR/mesa_black.bim
## Replace "---" chromosome numbers with something plink likes better (to be overwritten below)
#
#plink --bfile $DESTDIR/mesa_black \
#	--update-chr ../data/processed/snp_annotations/snp_annot_hg19_nodups.txt 1 3 \
#	--update-map ../data/processed/snp_annotations/snp_annot_hg19_nodups.txt 2 3 \
#	--maf 0.001 \
#	--make-bed \
#	--out $DESTDIR/mesa_black
#
#plink --bfile $DESTDIR/mesa_black \
#	--chr 1-22 \
#	--make-bed \
#	--out $DESTDIR/mesa_black
