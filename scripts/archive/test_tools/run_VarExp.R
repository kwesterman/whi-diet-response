library(tidyverse)
library(VarExp)
library(GWASTools)


bim_to_SnADF <- function(bfile) {
  bim_df <- read_tsv(paste0(bfile, ".bim"), 
                     col_names=c("chromosome", "rsID", "cm", 
                                 "position", "A1", "A2")) %>%
    mutate(snpID=1:nrow(.)) %>%
    mutate_at(vars(chromosome, position), as.integer)
  SnpAnnotationDataFrame(data.frame(bim_df, stringsAsFactors=F))
}

fam_to_ScADF <- function(bfile) {
  phenos <- read_delim("../data/processed/gen3/whi_white_gwas_phenos.txt", 
                       delim=" ") %>%
    mutate(ldl=ifelse(lipid_med, ldl / 0.75, ldl),
           glu=ifelse(dm_med, glu / 0.75, glu),
           sbp=ifelse(ht_med, sbp + 15, sbp))
  fam_df <- read_delim(paste0(bfile, ".fam"), 
                       delim=" ", col_names=c("FID", "IID", "father", "mother", 
                                              "sex", "pheno")) %>%
    left_join(phenos, by=c("IID")) %>%
    mutate(scanID=IID)
  ScanAnnotationDataFrame(data.frame(fam_df, stringsAsFactors=F))
}

bfile <- "../data/processed/whi_subsets/whi_bmi_suggestive"
snp_annot <- getAnnotation(bim_to_SnADF(bfile))
ss_annot <- read_tsv("../data/processed/gen3/sfa_bmi/sfa_bmi_suggestive.res_raw") %>%
  inner_join(snp_annot, by="snpID")

GWAS <- ss_annot %>%
  select(RSID=rsID, CHR=chromosome, POS=position, A0=A1, FREQ_A0=MAF,
         MAIN_EFFECT=Est, INT_EFFECT=GxE.Est)

phenos <- getAnnotation(fam_to_ScADF(bfile))
COHORT <- data.frame(
  Cohort=1,
  PHENO_N=sum(!is.na(phenos$bmi)),
  PHENO_Mean=mean(phenos$bmi, na.rm=T),
  PHENO_SD=sd(phenos$bmi, na.rm=T),
  EXPO_Mean=mean(phenos$sfa, na.rm=T),
  EXPO_SD=sd(phenos$sfa, na.rm=T)
)


C <- getGenoCorMatrix(GWAS$RSID, GWAS$CHR, GWAS$POS, GWAS$A0, "EUR",
                      pruning=T, rthresh=0.4)
GWAS <- checkInput(GWAS, colnames(C))

parsY <- calculateParamsFromIndParams(COHORT$PHENO_N, COHORT$PHENO_Mean, COHORT$PHENO_SD)
parsE <- calculateParamsFromIndParams(COHORT$PHENO_N, COHORT$EXPO_Mean, COHORT$EXPO_SD)

std_betaG <- standardizeBeta(GWAS$MAIN_EFFECT, GWAS$INT_EFFECT, GWAS$FREQ_A0, parsE[1], parsE[2], type = "G")
std_betaI <- standardizeBeta(GWAS$MAIN_EFFECT, GWAS$INT_EFFECT, GWAS$FREQ_A0, parsE[1], parsE[2], type = "I")

fracG <- calculateVarExp(std_betaG, std_betaI, C, parsY[2], sum(COHORT$PHENO_N), "G")
fracI <- calculateVarExp(std_betaG, std_betaI, C, parsY[2], sum(COHORT$PHENO_N), "I")
fracJ <- calculateVarExp(std_betaG, std_betaI, C, parsY[2], sum(COHORT$PHENO_N), "J")

print(fracG)
print(fracI)
print(fracJ)

# annotate_sumstats <- function(ss, bfile, maf_filter=0.01) {
#   anno <- bim_to_SnADF(bfile) %>%
#     getAnnotation() %>%
#     select(snpID, rsID, chromosome, position, A1, A2)
#   ss %>%
#     filter(MAF>maf_filter) %>%
#     select(snpID, MAF, n, GxE.Est, GxE.SE, GxE.pval) %>%
#     inner_join(anno, by="snpID") %>%
#     select(SNP=rsID, CHR=chromosome, BP=position, A1, A2, MAF, N=n, 
#            BETA=GxE.Est, SE=GxE.SE, P=GxE.pval)
# }




