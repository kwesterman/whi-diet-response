#!/bin/bash

module load plink
module load plink2

DATADIR=/cluster/home/kweste01/kw/ncbi/mesa/MESA_SHARe_imputation_CHN_1kg_10032012
DESTDIR=/cluster/home/kweste01/kw/diet_response/int/mesa

mkdir -p $DESTDIR

cat <(echo -e "ID_1 ID_2 missing\n0 0 0") \
	<(awk '{print $1,$1,0}' $DATADIR/imputation_sample_file_gen3_CHN.txt) \
	> $DESTDIR/mesa_chn.sample
cat $DATADIR/imputation_data/imp2_CHN_chr*_ref1kg > $DESTDIR/mesa_chn.gen

cat $DATADIR/SNP_info/imp2_CHN_chr*_ref1kg_info > $DESTDIR/mesa_chn.snpinfo
awk '$5>0.3 {print $2}' $DESTDIR/mesa_chn.snpinfo > $DESTDIR/mesa_chn_good_snps.txt

awk -F, '{print $1,$1}' ../int/metaData_mesa.csv > $DESTDIR/mesa_ids.txt

plink2 --gen $DESTDIR/mesa_chn.gen \
	--allow-extra-chr \
	--sample $DESTDIR/mesa_chn.sample \
	--extract $DESTDIR/mesa_chn_good_snps.txt \
	--keep $DESTDIR/mesa_ids.txt \
	--sort-vars \
	--make-pgen \
	--out $DESTDIR/mesa_chn

plink2 --pfile $DESTDIR/mesa_chn \
	--allow-extra-chr \
	--make-bed \
	--out $DESTDIR/mesa_chn

plink --bfile $DESTDIR/mesa_chn \
	--allow-extra-chr \
	--update-chr ../int/snp_annotations/snp_annot_hg19_nodups.txt 1 3 \
	--update-map ../int/snp_annotations/snp_annot_hg19_nodups.txt 2 3 \
	--maf 0.001 \
	--make-bed \
	--out $DESTDIR/mesa_chn
