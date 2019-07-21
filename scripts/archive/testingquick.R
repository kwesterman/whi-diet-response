ldl_phen <- process_metadata(clean_phenos_whi, "myhei", "ldl", "white")
ldl_wg <- add_genotypes(ldl_phen, "../data/processed/whi_ldl_suggestive.rds")
ldl_res <- cs_trial_compare(ldl_wg, "sfa", "ldl")

sbp_phen <- process_metadata(clean_phenos_whi, "myhei", "sbp", "white")
sbp_wg <- add_genotypes(sbp_phen, "../data/processed/whi_bp_suggestive.rds")
sbp_res <- cs_trial_compare(sbp_wg, "myhei", "sbp")

bmi_phen <- process_metadata(clean_phenos_whi, "myhei", "bmi", "white")
bmi_wg <- add_genotypes(bmi_phen, "../data/processed/whi_bmi_genomewide16.rds")
bmi_res <- cs_trial_compare(bmi_wg, "myhei", "bmi")


sbp_res_rfadj <- cs_trial_compare(sbp_wg, "myhei", "sbp")
sbp_res_doubleadj <- cs_trial_compare(sbp_wg, "myhei", "sbp")
sbp_res_prodadj <- cs_trial_compare(sbp_wg, "myhei", "sbp")


my_compare <- function(df1, df2, thresh=1) {
  compare_df <- data.frame(a=df1$estimate, b=df2$estimate) %>%
    filter(df1$p.value < thresh)
  plot(compare_df$a, compare_df$b)
  cor.test(~ a + b, data=compare_df)
}


n3_phen <- process_metadata(clean_phenos_whi, "n3", "tg", "white")
n3_wg <- add_genotypes(n3_phen, "../data/processed/whi_n3tg_response.rds")
n3_res <- cs_trial_compare(n3_wg, "n3", "tg")




ukbb <- read_table("../data/processed/gen2/sfa_sbp/ukbb_res.txt.gz")
ukbb <- mutate(ukbb, snp=paste0(gsub("_.*", "", SNP), "_", A1))




sbp_phen_nomeds <- process_metadata(filter(clean_phenos_whi, ht_med == F), 
                                    "myhei", "sbp", "white")
sbp_wg_nomeds <- add_genotypes(sbp_phen_nomeds, "../data/processed/whi_bp_suggestive.rds")
sbp_res_nomeds <- cs_trial_compare(sbp_wg_nomeds, "myhei", "sbp")




ldl_wg_nostatin <- list(long=ldl_wg$long, baseline=filter(ldl_wg$baseline, lipid_med == F))
ldl_res_nostatin <- cs_trial_compare(ldl_wg_nostatin, "sfa", "ldl")




