#!/bin/bash

PHENO=$1
module load metal

echo """
MARKER SNP 
ALLELE A1 A2
EFFECT BETA
PVALUE P
WEIGHT NMISS
SCHEME SAMPLESIZE
PROCESS ../data/processed/gen6/fat_${PHENO}_whi_white.res_annot

MARKER SNP 
ALLELE A1 A2
EFFECT BETA
PVALUE P
WEIGHT NMISS
PROCESS ../data/processed/ukbb_res/${PHENO}/fat_${PHENO}.res

OUTFILE ${PHENO}_MA_ .tbl
MINWEIGHT 125000
ANALYZE

QUIT
""" | metal
