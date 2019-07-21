# a <- lapply(1:22, function(chrom) {
  # read_table2(paste0("../data/processed/ukbb_res/bmi/PLINK_all_assoc_INFOcut30_MAF001_bmiint_CHR",
  #                    chrom, ".assoc.linear"))
# })
# saveRDS(a, file="my_ukbb_save.rds")

bmi_ss_list <- list()

# for (chrom in 1:22) {
#   print(chrom)
#   bmi_ss_list[[chrom]] <- read_table2(
#     # paste0("../data/processed/ukbb_res/bmi/PLINK_all_assoc_INFOcut30_MAF001_bmiint_CHR",
#     #        chrom, ".assoc.linear"),
#     "../data/processed/ukbb_res/bmi/plink_res_concat.assoc.linear",
#     col_types=cols_only(CHR="i", SNP="c", BP="i", TEST="c", BETA="d", P="d"))
# }

bmi_gwis_res <- read_table2(
      "../data/processed/ukbb_res/bmi/plink_res_concat.assoc.linear",
      col_types=cols_only(CHR="i", SNP="c", BP="i", TEST="c", BETA="d", P="d"))

saveRDS(bmi_gwis_res, file="my_ukbb_save.rds")

bmi_gwis_res$SNP_with_alleles <- bmi_gwis_res$SNP
bmi_gwis_res$SNP <- gsub("_.*", "", bmi_gwis_res$SNP_with_alleles)

library(qqman)
manhattan(filter(bmi_gwis_res, P < 0.01), chr="CHR", bp="BP", snp="SNP")
