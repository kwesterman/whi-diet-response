library(tidyverse)

whi_event_names <- c("CHD", "MI", "STRKISCH", "STRKHEMO", "DEATH", "ANYCANCER")#, "CREVASC")
whi_event_days <- paste0(whi_event_names, "DY")
cvd_death_causes <- c(11:19)
outcome_data_c1 <- read_tsv("../data/raw/whi/phen/outcomes_ctos_c1.txt", skip=10) %>%
  select(SUBJID, whi_event_names, whi_event_days, DEATHCAUSE, ENDFOLLOWDY)
outcome_data_c2 <- read_tsv("../data/raw/whi/phen/outcomes_ctos_c2.txt", skip=10) %>%
  select(SUBJID, whi_event_names, whi_event_days, ENDFOLLOWDY)
outcome_data <- bind_rows(outcome_data_c1, outcome_data_c2) %>%
  mutate(chd=CHD,
         time_to_chd=pmin(CHDDY, ENDFOLLOWDY, na.rm=T),
         isch_stroke=STRKISCH,
         time_to_isch_stroke=pmin(STRKISCHDY, ENDFOLLOWDY, na.rm=T),
         hemo_stroke=STRKHEMO,
         time_to_hemo_stroke=pmin(STRKHEMODY, ENDFOLLOWDY, na.rm=T),
         death=DEATH,
         time_to_death=pmin(DEATHDY, ENDFOLLOWDY, na.rm=T),
         cvd_death=DEATH * (DEATHCAUSE %in% cvd_death_causes),
         time_to_cvd_death=pmin(DEATHDY, ENDFOLLOWDY, na.rm=T),
         non_cvd_death=DEATH * !(DEATHCAUSE %in% cvd_death_causes),
         time_to_non_cvd_death=pmin(DEATHDY, ENDFOLLOWDY, na.rm=T),
         cancer=ANYCANCER,
         time_to_cancer=pmin(ANYCANCERDY, ENDFOLLOWDY, na.rm=T),
         mi=MI,
         time_to_mi=pmin(MIDY, ENDFOLLOWDY, na.rm=T))

dm_data_c1 <- read_tsv("../data/raw/whi/phen/outcomes_self_c1.txt", skip=10)
dm_data_c2 <- read_tsv("../data/raw/whi/phen/outcomes_self_c2.txt", skip=10)
dm_data <- bind_rows(dm_data_c1, dm_data_c2) %>%
  mutate(dm=F33DIABPILLS | F33DIABINSULIN | F33DIABTX,
         time_to_dm=pmin(F33DIABPILLSDY, F33DIABINSULINDY, F33DIABTXDY, na.rm=T)) %>%
  select(SUBJID, dm, time_to_dm)

outcome_data <- outcome_data %>%
  left_join(dm_data, by="SUBJID") %>%
  mutate(time_to_dm=pmin(time_to_dm, ENDFOLLOWDY, na.rm=T)) %>%
  select(subjID=SUBJID, chd, time_to_chd, isch_stroke, time_to_isch_stroke, 
         hemo_stroke, time_to_hemo_stroke, death, time_to_death,
         cvd_death, time_to_cvd_death, non_cvd_death, time_to_non_cvd_death,
         dm, time_to_dm, mi, time_to_mi,
         cancer, time_to_cancer)

write_csv(outcome_data, "../data/processed/outcomes_whi.csv")
