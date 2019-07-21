library(tidyverse)
library(GWASTools)
library(SNPRelate)

bim_to_SnADF <- function(bfile) {
  bim_df <- read_tsv(paste0(bfile, ".bim"), 
                     col_names=c("chromosome", "rsID", "cm", 
                                 "position", "A1", "A2")) %>%
    mutate(snpID=1:nrow(.)) %>%
    mutate_at(vars(chromosome, position), as.integer)
  SnpAnnotationDataFrame(data.frame(bim_df, stringsAsFactors=F))
}

fam_to_ScADF <- function(bfile) {
  phenos <- read_delim("../data/processed/gen4/fhs_gwas_phenos.txt", 
                       delim=" ") %>%
    mutate(myvar=sfa>22)
    # mutate(myhei=-scale(palm) + scale(pufa) + scale(n3))
  fam_df <- read_delim(paste0(bfile, ".fam"), 
                       delim="\t", col_names=c("FID", "IID", "father", "mother", 
                                               "sex", "pheno")) %>%
    select(IID) %>%
    left_join(phenos, by=c("IID")) %>%
    mutate(scanID=IID)
  ScanAnnotationDataFrame(data.frame(fam_df, stringsAsFactors=F))
}

make_gds <- function(bfile, gds_name, summary=F) {
  snpgdsBED2GDS(paste0(bfile, ".bed"),
                paste0(bfile, ".fam"),
                paste0(bfile, ".bim"),
                gds_name,
                cvt.snpid="int")
  if (summary) snpgdsSummary(gds_name)
}


bfile <- paste0("../data/processed/fhs_subsets/fhs_bmi_hypothesis")
gds_name <- paste0(bfile, ".gds")
make_gds(bfile, gds_name)
gds <- GdsGenotypeReader(openfn.gds(gds_name, allow.fork=T))

snpAnnot <- bim_to_SnADF(bfile)
scanAnnot <- fam_to_ScADF(bfile)
genoData <- GenotypeData(gds, snpAnnot=snpAnnot, scanAnnot=scanAnnot)

sfa_only <- c("sfa", "age")
sfa_carbEx_covars <- c("sfa", "mufa", "pufa", "pro", "alc", "tot_cal", "age")
sfa_pufaEx_covars <- c("sfa", "mufa", "pro", "cho", "alc", "tot_cal", "age")
fat_carbEx_covars <- c("fat", "pro", "alc", "tot_cal", "age")
# covars <- c("myhei", "tot_cal", "age")
outcome <- "logBMI"
ivar <- "sfa"
phenos <- getAnnotation(scanAnnot)

run_test <- function(covars, ivar, robust) {
  a <- assocRegression(
    genoData,
    outcome=outcome,
    model.type="linear",
    gene.action="additive",
    covar=covars,
    ivar=ivar,
    robust=robust)
  a %>%
    select(snpID, GxE.Est, GxE.pval, JOINT=Joint.pval) #%>% 
    # mutate_at(vars(-snpID), round, 3)
}

run_test(c("myhei", "tot_cal", "age"), "myvar", robust=T)


lapply(list(sfa_only, sfa_carbEx_covars, sfa_pufaEx_covars),
       run_test, "sfa", robust=T)
