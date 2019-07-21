# For at least testing GWAS through 

library(GWASTools)
library(SNPRelate)
library(tidyverse)
library(foreach)
library(doParallel)

bfile <- "whi_white_nonDM" 

# SNP annotations
bim_df <- read_tsv(paste0("../data/processed/gen3/", bfile_prefix, ".bim"), 
                   col_names=c("chromosome", "rsID", "cm", 
                               "position", "A1", "A2")) %>%
  mutate(snpID=1:nrow(.)) %>%
  mutate_at(vars(chromosome, position), as.integer)
snpAnnot <- SnpAnnotationDataFrame(data.frame(bim_df, stringsAsFactors=F))

# Scan annotations
phenos <- read_delim("../data/processed/gen3/whi_white_gwas_phenos.txt", 
                     delim=" ")
fam_df <- read_delim(paste0("../data/processed/gen3/", bfile_prefix, ".fam"), 
                     delim=" ", col_names=c("FID", "IID", "father", "mother", 
                                            "sex", "pheno")) %>%
  left_join(phenos, by=c("IID")) %>%
  mutate(scanID=IID)
scanAnnot <- ScanAnnotationDataFrame(data.frame(fam_df, stringsAsFactors=F))

# Genotype data (hard calls for the moment)
# snpgdsBED2GDS(paste0("../data/processed/gen3/", bfile_prefix, ".bed"),
#               paste0("../data/processed/gen3/", bfile_prefix, ".fam"),
#               paste0("../data/processed/gen3/", bfile_prefix, ".bim"),
#               "current_gds",
#               cvt.snpid="int")
# snpgdsSummary("current_gds")
gds <- GdsGenotypeReader(openfn.gds("current_gds", allow.fork=T))
genoData <- GenotypeData(gds, snpAnnot=snpAnnot, scanAnnot=scanAnnot)

# GWAS
outcome <- "logGLU"
covars <- c("age", "fat", "pro", "alc", "tot_cal", paste0("PC", 1:5))
ivar <- "fat"

run_gwas_chunk <- function(outcome, covars, ivar, idx, robust=F) {
  assocRegression(
    genoData,
    outcome=outcome,
    model.type="linear",
    covar=covars,
    ivar=ivar,
    robust=robust,
    snpStart=min(which(buckets == idx)),
    snpEnd=max(which(buckets == idx)))
}

num_cores <- 16
buckets <- as.integer(cut(1:nsnp(genoData), num_cores))
cl <- makeForkCluster(num_cores)
registerDoParallel(cl)
all_res <- foreach(idx=unique(buckets), .combine=rbind, .packages="GWASTools") %dopar%
  run_gwas_chunk(idx, robust=T)
stopCluster(cl)

saveRDS(all_res, "../data/processed/gen3/whi_robust_fat_glu.rds")

# all_results <- lapply(unique(buckets), function(grp) {
#   assocRegression(
#     genoData,
#     outcome="logGLU",
#     model.type="linear",
#     covar=c("age", "fat", "pro", "alc", 
#             "tot_cal", paste0("PC", 1:5)),
#     ivar="fat",
#     robust=F,
#     snpStart=min(which(buckets == grp)),
#     snpEnd=max(which(buckets == grp)))
# })
# 
# test_results <- assocRegression(
#   genoData,
#   outcome="logGLU",
#   model.type="linear",
#   covar=c("age", "fat", "pro", "alc", 
#           "tot_cal", paste0("PC", 1:5)),
#   ivar="fat",
#   robust=F,
#   snpStart=1,
#   snpEnd=1000)
# 
# test_results <- assocRegression(
#   genoData,
#   outcome="logGLU",
#   model.type="linear",
#   covar=c("age", "fat", "pro", "alc", 
#           "tot_cal", paste0("PC", 1:5)),
#   ivar="fat",
#   PPLtest=F,
#   robust=F,
#   snpStart=1,
#   snpEnd=1000)

# system.time(test_results_robust <- assocRegression(
#   genoData,
#   outcome="logGLU",
#   model.type="linear",
#   covar=c("age", "fat", "pro", "alc", 
#           "tot_cal", paste0("PC", 1:5)),
#   ivar="fat",
#   robust=T))

make_qqplot <- function(p_vec, plotTitle="Title") {
  p_vec <- p_vec[!is.na(p_vec)]
  qqplot(-log10(1:length(p_vec) / length(p_vec)), -log10(p_vec), pch=".", 
         main=plotTitle, xlab="Expected (-logP)", ylab="Observed (-logP)")
  abline(0, 1, col="red")
}

gControl <- function(pVals) {
  # See van Iterson 2017 methods and/or Lehne 2015 code for details on genomic control for EWAS
  # Below is modeled after Lehne 2015
  lambda <- median(qchisq(pVals, df=1, lower.tail=F), na.rm=T)/qchisq(0.5, df=1)
  round(lambda, 2)
}

print(gControl(all_res$GxE.pval))

# print("Non-robust:")
# print(gControl(test_results$GxE.pval))
# print("Robust:")
# print(gControl(test_results_robust$GxE.pval))

