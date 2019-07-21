#!/bin/bash

module load plink2

DATADIR=/cluster/home/kweste01/kw/diet_response/data/raw/bprhs/gen
DESTDIR=/cluster/home/kweste01/kw/diet_response/data/processed/bprhs
D2P=/cluster/home/kweste01/kw/opt/dose2plink

# Loop through chromosomes, using dose2plink to transpose .dose files
#for ((i=1; i<=22; i++)); do 
#	awk '{print $3,$4,$5,$6,$7,$8,$9}' <(tail -n +2 $DATADIR/CHR${i}.ANNOT.txt) > $DESTDIR/tmpinfobody_chr${i}.txt
#	cat <(echo "SNP Al1 Al2 Freq1 MAF Qual Rsq") $DESTDIR/tmpinfobody_chr${i}.txt > $DESTDIR/tmpinfo_chr${i}.txt
#	awk '{print $1"->"$1, "DOSE"}' $DATADIR/chr${i}.dose > $DESTDIR/tmpdoseidcols_chr${i}.txt
#	paste $DESTDIR/tmpdoseidcols_chr${i}.txt <(cut -f 2- $DATADIR/chr${i}.dose) | tail -n +2 > $DESTDIR/tmpdose_chr${i}.dose
#	$D2P -dose $DESTDIR/tmpdose_chr${i}.dose -info $DESTDIR/tmpinfo_chr${i}.txt -out $DESTDIR/bprhs_chr${i}
#done

# Loop back through chromosomes, tracking low-quality SNPs and unzipping .pdat files
#gunzip < $DESTDIR/bprhs_chr1.pdat.gz | head -1 > $DESTDIR/bprhs.pdat
#for ((i=1; i<=22; i++)); do
#	awk '{if ($9<0.3) print $3}' $DATADIR/CHR${i}.ANNOT.txt > $DESTDIR/chr${i}_low_qual_snps.txt
#	gunzip < $DESTDIR/bprhs_chr${i}.pdat.gz | tail -n +2 >> $DESTDIR/bprhs.pdat
#done
#
#cat $DESTDIR/chr*_low_qual_snps.txt > $DESTDIR/low_qual_snps.txt

SNPANNO=../data/processed/snp_annotations/snp_annot_hg19_nodups.txt  # Contains 1000 Genomes SNP locations, rsIDs, and alleles from Ensembl
plink2 --import-dosage $DESTDIR/bprhs.pdat \
	--psam $DESTDIR/bprhs_chr1.pfam \
	--exclude $DESTDIR/low_qual_snps.txt \
	--ref-allele $SNPANNO 4 3 \
	--maf 0.001 \
	--make-pgen \
	--out $DESTDIR/bprhs
