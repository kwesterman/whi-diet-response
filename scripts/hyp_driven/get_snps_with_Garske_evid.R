library(tidyverse)
library(readxl)
library(GenomicRanges)


# Lipid-responsive promoter or enhancer regions in 
garske_promoter_loci <- read_excel("../data/raw/literature/lipid_responsive_promoter_peaks.xlsx")
garske_enhancer_loci <- read_excel("../data/raw/literature/lipid_responsive_enhancer_peaks.xlsx")
garske_loci <- bind_rows(garske_promoter_loci, garske_enhancer_loci)
garske_ranges <- GRanges(seqnames=garske_loci$peakChr, 
                         ranges=IRanges(start=garske_loci$peakStart, 
                                        end=garske_loci$peakEnd))

# SNPs tested in BMI GWIS
bmi_snp_loci <- read_delim("../data/processed/whi_prediction/bmi_tested_variant_loci.txt",
                           delim=" ") %>%
  filter(!is.na(bp))
all_snp_ranges <- GRanges(seqnames=bmi_snp_loci$chr,
                          ranges=IRanges(start=bmi_snp_loci$bp,
                                         end=bmi_snp_loci$bp))

# Find overlaps to prioritize SNPs
overlaps <- findOverlaps(query=garske_ranges, subject=all_snp_ranges)
bmi_priority_snps <- bmi_snp_loci$SNP[subjectHits(overlaps)]
write(bmi_priority_snps, "../data/processed/whi_prediction/bmi_prioritized_variants.txt")
