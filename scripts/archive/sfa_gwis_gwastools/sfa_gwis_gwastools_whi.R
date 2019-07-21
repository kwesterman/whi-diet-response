silent <- lapply(
c("knitr", "tidyverse", "cowplot", "doParallel", "kableExtra",
"glmnet", "broom", "GWASTools", "SNPRelate"), library, character.only=T)


args <- commandArgs(trailingOnly=T)
dv_withEx <- args[1]
rf <- args[2]
main_effect_threshold <- args[3]

INT <- function(x) qnorm((rank(x, na.last="keep") - 0.5) / sum(!is.na(x)))

winsorize <- function(x, num_SDs=5) {
  bounds <- mean(x, na.rm=T) + num_SDs * c(-1, 1) * sd(x, na.rm=T)
  case_when(x < bounds[1] ~ bounds[1],
            x > bounds[2] ~ bounds[2],
            TRUE ~ x)
}

make_qqplot <- function(p_vec, plotTitle="Title") {
  p_vec <- p_vec[!is.na(p_vec)]
  qqplot(-log10(1:length(p_vec) / length(p_vec)), -log10(p_vec), pch=".", 
         main=plotTitle, xlab="Expected (-logP)", ylab="Observed (-logP)")
  abline(0, 1, col="red")
}

gControl <- function(pVals) {
  # See van Iterson 2017 methods and/or Lehne 2015 code for details on genomic control for EWAS
  # Below is modeled after Lehne 2015
  lambda <- median(qchisq(pVals, df=1, lower.tail=F), na.rm=T) / qchisq(0.5, df=1)
  round(lambda, 2)
}

bim_to_SnADF <- function(bfile) {
  bim_df <- read_tsv(paste0(bfile, ".bim"), 
                     col_names=c("chromosome", "rsID", "cm", 
                                 "position", "A1", "A2")) %>%
    mutate(snpID=1:nrow(.)) %>%
    mutate_at(vars(chromosome, position), as.integer)
  SnpAnnotationDataFrame(data.frame(bim_df, stringsAsFactors=F))
}

fam_to_ScADF <- function(bfile) {
  phenos <- read_delim("../data/processed/gen4/whi_white_gwas_phenos.txt", 
                       delim=" ")
  fam_df <- read_delim(paste0(bfile, ".fam"), 
                       delim=" ", col_names=c("FID", "IID", "father", "mother", 
                                              "sex", "pheno")) %>%
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

run_gwis_chunk <- function(genoData, outcome, covars, ivar, start, stop, robust=F) {
  assocRegression(
    genoData,
    outcome=outcome,
    model.type="linear",
    covar=covars,
    ivar=ivar,
    robust=robust,
    snpStart=start,
    snpEnd=stop)
}

covar_sets <- list(
  fat_carbEx=c("age", "fat", "pro", "alc", "tot_cal", paste0("PC", 1:5)),
  sfa_carbEx=c("age", "sfa", "mufa", "pufa", "pro", "alc", "tot_cal", 
               paste0("PC", 1:5)),
  myhei=c("age", "myhei", paste0("PC", 1:5)),
  fat_carbEx_Bin=c("age", "fatBin", "pro", "alc", "tot_cal", paste0("PC", 1:5)),
  sfa_carbEx_Bin=c("age", "sfaBin", "mufa", "pufa", "pro", "alc", "tot_cal", 
               paste0("PC", 1:5)),
  myhei_Bin=c("age", "myheiBin", paste0("PC", 1:5)),
  sfa_binary=c("sfa_binary", "age", "tot_cal", paste0("PC", 1:5))
)

interaction_vars <- list(
  fat_carbEx="fat",
  sfa_carbEx="sfa",
  myhei="myhei",
  fat_carbEx_Bin="fatBin",
  sfa_carbEx_Bin="sfaBin",
  myhei_Bin="myheiBin",
  sfa_binary="sfa_binary"
)

outcome_transforms <- list(
  bmi="logBMI",
  # hsCRP="logHSCRP",
  tg="logTG",
  glu="logGLU",
  ldl="ldl",
  hdl="logHDL",
  sbp="logSBP",
  delta_bmi="delta_bmi",
  delta_tg="delta_tg",
  delta_glu="delta_glu",
  delta_ldl="delta_ldl",
  delta_sbp="delta_sbp"
)

outcome_basenames <- list(
  bmi="bmi",
  # hsCRP="logHSCRP",
  tg="tg",
  glu="glu",
  ldl="ldl",
  hdl="hdl",
  sbp="sbp",
  delta_bmi="bmi",
  delta_tg="tg",
  delta_glu="glu",
  delta_ldl="ldl",
  delta_sbp="sbp"
)

run_gwis <- function(outcome, dv, genoset, robust=T) {
    covars <- covar_sets[[dv]]
    ivar <- interaction_vars[[dv]]
    outcome_transform <- outcome_transforms[[outcome]]
    outcome_basename <- outcome_basenames[[outcome]]
    bfile <- paste0("../data/processed/whi_subsets/whi_", 
                    "bmi",  # TEMPORARY HACK
                    "_", genoset)
    snpAnnot <- bim_to_SnADF(bfile)
    scanAnnot <- fam_to_ScADF(bfile)
    gds_name <- paste0(bfile, ".gds")
    make_gds(bfile, gds_name)
    gds <- GdsGenotypeReader(openfn.gds(gds_name, allow.fork=T))
    genoData <- GenotypeData(gds, snpAnnot=snpAnnot, scanAnnot=scanAnnot)
    
    num_cores <- detectCores()
    chunks <- as.integer(cut(1:nsnp(genoData), num_cores))
    cl <- makeForkCluster(num_cores)
    registerDoParallel(cl)
    res <- foreach(idx=unique(chunks), .combine=rbind, 
                       .packages="GWASTools") %dopar%
      run_gwis_chunk(genoData, outcome_transform, covars, ivar, 
                     min(which(chunks == idx)), max(which(chunks == idx)), 
                     robust=robust)
    stopCluster(cl)
    
    close(gds)
    
    res
}

####### RUN GWIS #######

gwis_res <- run_gwis(rf, dv_withEx, main_effect_threshold)

annotate_sumstats <- function(ss, bfile, maf_filter=0.01) {
  anno <- bim_to_SnADF(bfile) %>%
    getAnnotation() %>%
    select(snpID, rsID, chromosome, position, A1, A2)
  ss %>%
    filter(MAF > maf_filter) %>%
    select(snpID, MAF, n, GxE.Est, GxE.SE, GxE.pval) %>%
    inner_join(anno, by="snpID") %>%
    select(SNP=rsID, CHR=chromosome, BP=position, A1, A2, MAF, N=n,
           BETA=GxE.Est, SE=GxE.SE, P=GxE.pval)
}

dv <- interaction_vars[[dv_withEx]]

write_tsv(gwis_res, 
          paste0("../data/processed/gen4/whi_res/whi_white_", dv, "_", rf, "_",
                 main_effect_threshold, ".res_raw"))
res_anno <- annotate_sumstats(gwis_res,
                              paste0("../data/processed/whi_subsets/whi_",
                                     rf, "_", main_effect_threshold))
write_tsv(res_anno, 
          paste0("../data/processed/gen4/whi_res/whi_white_", dv, "_", rf, "_",
                 main_effect_threshold, ".res"))