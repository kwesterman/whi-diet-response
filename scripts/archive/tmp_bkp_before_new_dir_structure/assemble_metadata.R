# Assemble metadata (demographics, CVD events, lab values, etc.) on subjects

suppressMessages(silent <- lapply(c("tidyverse","foreign"), 
                                  library, character.only=T))

# FHS --------------------------------------------------------------------------

# Risk factors

phenos_fhs_c1 <- read_tsv("../data/fhs/phen/phenos_fhs_c1.txt", skip=10)
phenos_fhs_c2 <- read_tsv("../data/fhs/phen/phenos_fhs_c2.txt", skip=10)
phenos_fhs <- bind_rows(phenos_fhs_c1, phenos_fhs_c2) %>%
  rename(subjID=shareid, sex=SEX, 
         age_5=AGE5, bmi_5=BMI5, 
         smk_now_5=CURRSMK5, cig_per_day_5=CPD5,  
         sbp_5=SBP5, glu_5=FASTING_BG5, chol_5=TC5, ldl_5=CALC_LDL5, hdl_5=HDL5, 
         tg_5=TRIG5, ht_med_5=HRX5, lipid_med_5=LIPRX5, dm_med_5=DMRX5,
         age_8=AGE8, bmi_8=BMI8, 
         smk_now_8=CURRSMK8, cig_per_day_8=CPD8,  
         sbp_8=SBP8, glu_8=FASTING_BG8, chol_8=TC8, ldl_8=CALC_LDL8, hdl_8=HDL8, 
         tg_8=TRIG8, ht_med_8=HRX8, lipid_med_8=LIPRX8, dm_med_8=DMRX8) %>%
  mutate(sex=ifelse(sex==1, "M", "F"),
         race="white") %>%
  mutate_at(vars(ht_med_5, lipid_med_5, dm_med_5, 
                 ht_med_8, lipid_med_8, dm_med_8), as.logical) %>%
  select(subjID, sex, race,
         age_5, bmi_5, smk_now_5, cig_per_day_5, 
         sbp_5, glu_5, chol_5, ldl_5, hdl_5, tg_5, 
         ht_med_5, lipid_med_5, dm_med_5,
         age_8, bmi_8, smk_now_8, cig_per_day_8, 
         sbp_8, glu_8, chol_8, ldl_8, hdl_8, tg_8, 
         ht_med_8, lipid_med_8, dm_med_8)

crp_fhs_c1 <- read_tsv("../data/fhs/phen/crp_c1.txt", skip=10)
crp_fhs_c2 <- read_tsv("../data/fhs/phen/crp_c2.txt", skip=10)
crp_fhs <- bind_rows(crp_fhs_c1, crp_fhs_c2) %>%
  rename(subjID=shareid, hscrp_8=crp)

questionnaire_fhs_c1 <- read_tsv("../data/fhs/phen/exam8Data_fhs_c1.txt", skip=10)
questionnaire_fhs_c2 <- read_tsv("../data/fhs/phen/exam8Data_fhs_c2.txt", skip=10)
questionnaire_fhs <- bind_rows(questionnaire_fhs_c1, questionnaire_fhs_c2) %>%
  rename(subjID=shareid, cig_start_age=H065) %>%
  select(subjID, cig_start_age)

phenoData_fhs <- phenos_fhs %>%
  left_join(crp_fhs, by="subjID") %>%
  left_join(questionnaire_fhs, by="subjID") %>%
  mutate(subjID=as.character(subjID),
         smk_py_5=ifelse(smk_now_5,  # Very rough -- assumes no one has quit
                         cig_per_day_5 / 20 * (age_5 - cig_start_age), 0),
         smk_py_8=ifelse(smk_now_8, 
                         cig_per_day_8 /20 * (age_8 - cig_start_age), 0))

# Diet

ffq_fhs_ex5_c1 <- read_tsv("../data/fhs/diet/ffq_ex5_c1.txt", skip=10)
ffq_fhs_ex5_c2 <- read_tsv("../data/fhs/diet/ffq_ex5_c2.txt", skip=10)
ffq_fhs_ex5 <- bind_rows(ffq_fhs_ex5_c1, ffq_fhs_ex5_c2) %>%
  mutate(subjID=as.character(shareid),
         tot_cal_5=NUT_CALOR,
         sfa_5=NUT_SATFAT,
         mufa_5=NUT_MONFAT,
         pufa_5=NUT_POLY,
         animal_fat_5=NUT_AFAT,
         veg_fat_5=NUT_VFAT,
         carb_5=NUT_CARBO,
         sucrose_5=NUT_SUCR,
         sodium_5=NUT_SODIUM,
         fructose_5=NUT_FRUCT,
         n3_5=NUT_N3,
         n6_5=NUT_N6,
         linoleic_5=NUT_F182,
         palmitic_5=NUT_F160,
         nuts_5=FFD134) %>%
  mutate_at(vars(contains("_5"), -tot_cal_5), funs(. / tot_cal_5)) %>%
  select(subjID, contains("_5")) %>%
  mutate(sfa2pufa_5=sfa_5/pufa_5,
         n62n3_5=n6_5/n3_5)
ffq_fhs_ex8_c1 <- read_tsv("../data/fhs/diet/ffq_ex8_c1.txt", skip=10)
ffq_fhs_ex8_c2 <- read_tsv("../data/fhs/diet/ffq_ex8_c2.txt", skip=10)
ffq_fhs_ex8 <- bind_rows(ffq_fhs_ex8_c1, ffq_fhs_ex8_c2) %>%
  mutate(subjID=as.character(shareid),
         tot_cal_8=NUT_CALOR,
         sfa_8=NUT_SATFAT,
         mufa_8=NUT_MONFAT,
         pufa_8=NUT_POLY,
         animal_fat_8=NUT_AFAT,
         veg_fat_8=NUT_VFAT,
         carb_8=NUT_CARBO,
         sucrose_8=NUT_SUCR,
         sodium_8=NUT_SODIUM,
         fructose_8=NUT_FRUCT,
         n3_8=NUT_PFN302,
         n6_8=NUT_N602,
         linoleic_8=NUT_F182,
         palmitic_8=NUT_F160,
         nuts_8=FFD134) %>%
  mutate_at(vars(contains("_8"), -tot_cal_8), funs(. / tot_cal_8)) %>%
  select(subjID, contains("_8")) %>%
  mutate(sfa2pufa_8=sfa_8/pufa_8,
         n62n3_8=n6_8/n3_8)
ffq_fhs <- full_join(ffq_fhs_ex5, ffq_fhs_ex8, by="subjID")

fhsMetaData <- inner_join(phenoData_fhs, ffq_fhs, by="subjID")
saveRDS(fhsMetaData, "../int/metaData_fhs.rds")
write_csv(fhsMetaData, "../int/metaData_fhs.csv")

# exam_dates_fhs_c1 <- read_tsv("../data/fhs/phen/exam_dates_fhs_c1.txt", skip=10, 
#                               col_types=cols_only(shareid="i",date8="i",date9="i"))
# exam_dates_fhs_c2 <- read_tsv("../data/fhs/phen/exam_dates_fhs_c2.txt", skip=10,
#                               col_types=cols_only(shareid="i",date8="i",date9="i"))
# exam_dates_fhs <- bind_rows(exam_dates_fhs_c1, exam_dates_fhs_c2)
# 
# soe2015_fhs_c1 <- read_tsv("../data/fhs/phen/soe2015_c1.txt", skip=10, col_types=cols(shareid="i"))
# soe2015_fhs_c2 <- read_tsv("../data/fhs/phen/soe2015_c2.txt", skip=10, col_types=cols(shareid="i"))
# soe2015_fhs <- bind_rows(soe2015_fhs_c1, soe2015_fhs_c2)
# 
# soe2015_fhs_clean <- inner_join(soe2015_fhs, exam_dates_fhs, by="shareid") %>%
#   filter(!is.na(date8),  # Had a visit during Exam 8 
#          EVENT %in% c(1:29)) %>%  # CVD events (my definition) and all death
#   mutate(eventType=case_when(EVENT %in% c(1:9,21:24) ~ "chd",  # MI variants, AP, CHD deaths
#                              EVENT %in% c(10:19,25) ~ "stroke",  # CVA variants, ABI, TIA, embolism, hemorrhage
#                              EVENT == 26 ~ "death_otherCVD",
#                              EVENT %in% 27:29 ~ "death_nonCVD"),
#          cvd=eventType %in% c("chd","stroke","death_otherCVD"),
#          time=DATE-date8) %>%
#   group_by(shareid) %>%
#   summarise(pastEvent=any(cvd==T & time<=0),  # Note if subject had an event before Exam 8
#             event=any(cvd==T & time>0),  # Future event if occurred after Exam 8
#             timeToEvent=ifelse(any(cvd==T & time>0), min(time[cvd==T & time>0]), NA),  # Earliest post-Exam 8 event time
#             eventType=ifelse(any(cvd==T & time>0), eventType[which.min(time[cvd==T & time>0])], as.character(NA)),
#             death=any(EVENT %in% 21:29),  # Did the person die?
#             timeToDeath=ifelse(any(death), time[EVENT %in% 21:29], NA),  # If they died, get time to death
#             incCHD=any(eventType=="chd" & time>0),
#             incStroke=any(eventType=="stroke" & time>0),
#             timeToExam9=median(date9-date8))  # Carry through for censorship times
# 
# surv2014_fhs_c1 <- read_tsv("../data/fhs/phen/survcvd2014_fhs_c1.txt", skip=10)
# surv2014_fhs_c2 <- read_tsv("../data/fhs/phen/survcvd2014_fhs_c2.txt", skip=10)
# surv2014_fhs <- bind_rows(surv2014_fhs_c1, surv2014_fhs_c2)
# 
# outcomeData_fhs <- left_join(surv2014_fhs, soe2015_fhs_clean, by="shareid") %>%
#   filter(shareid %in% sampleData_fhs$subjID) %>%
#   replace_na(list(event=F, pastEvent=F, death=F)) %>%  # Add "false" for events/deaths after left join
#   mutate(time=case_when(event==T ~ timeToEvent,  # Experienced an event after Exam 8 -> use that follow-up time
#                         death==T ~ timeToDeath,  # Didn't experience an event but died -> censor at death
#                         cvd==0 ~ cvddate,  # No event and no CVD in surv file -> censoring time from surv file
#                         TRUE ~ as.integer(timeToExam9))) %>%  # No event and CVD in surv file -> use exam 9 date for censor time
#   filter(!is.na(time)) %>%  # Remove those individuals with no Exam 9 date and thus no known censorship time
#   mutate(subjID=as.character(shareid)) %>%
#   select(subjID, pastEvent, event, eventType, time, incCHD, incStroke)
# ## NOTE: THE APPROACH ABOVE LEAVES A NUMBER OF SUBJECTS (114) WHO HAD A PAST EVENT, NO FUTURE EVENT,
# ## AND NO AVAILABLE EXAM 9 DATE WITH MISSING "TIME" VALUES AND THEY ARE THUS EXCLUDED
# ## THIS IS POTENTIALLY REASONABLE HERE BECAUSE NONE REPRESENT CASES (THE MORE CRITICAL SAMPLES TO KEEP)



# WHI --------------------------------------------------------------------------

# Risk factors

basicData_whi_c1 <- read_tsv("../data/whi/phen/basic_whi_c1.txt", skip=10)
basicData_whi_c2 <- read_tsv("../data/whi/phen/basic_whi_c2.txt", skip=10)
basicData_whi <- bind_rows(basicData_whi_c1, basicData_whi_c2) %>%
  rename(subjID=SUBJID, age=AGE, race=RACE) %>%
  mutate(sex="F",
         race=c("1"="amind", "2"="asian", "3"="black", "4"="hispanic", 
                "5"="white", "8"="other")[race]) %>%
  select(subjID, sex, age, race)

randomization_whi_c1 <- read_tsv("../data/whi/randomization_c1.txt", skip=10)
randomization_whi_c2 <- read_tsv("../data/whi/randomization_c2.txt", skip=10)
randomization_whi <- bind_rows(randomization_whi_c1, randomization_whi_c2) %>%
  rename(subjID=SUBJID) %>%
  mutate(dm_trial=as.logical(DMFLAG),
         dm_intervention=(DMARM == 1)) %>%
  select(subjID, dm_trial, dm_intervention)

behaviorData_whi_c1 <- read_tsv("../data/whi/phen/behavior_c1.txt", skip=10)
behaviorData_whi_c2 <- read_tsv("../data/whi/phen/behavior_c2.txt", skip=10)
behaviorData_whi <- bind_rows(behaviorData_whi_c1, behaviorData_whi_c2) %>%
  rename(subjID=SUBJID, smk_now=SMOKNOW, smk_py=PACKYRS) %>%
  select(subjID, smk_now, smk_py)

examData_whi_c1 <- read_tsv("../data/whi/phen/physical_whi_c1.txt", skip=10)
examData_whi_c2 <- read_tsv("../data/whi/phen/physical_whi_c2.txt", skip=10)
examData_whi <- bind_rows(examData_whi_c1, examData_whi_c2) %>%
  rename(subjID=SUBJID, visitYear=F80VY, sbp=SYST, bmi=BMIX) %>%
  select(subjID, visitYear, sbp, bmi)

labData1_whi_c1 <- read_tsv("../data/whi/phen/labs_c1.txt", skip=10)
labData1_whi_c2 <- read_tsv("../data/whi/phen/labs_c2.txt", skip=10)
labData1_whi <- bind_rows(labData1_whi_c1, labData1_whi_c2) %>%
  # filter(COREVTYP==1) %>%  # Only care about blood samples from first year
  filter(COREVTYP!=4) %>%  # Throw out non-routine (non-annual) visits
  rename(subjID=SUBJID, visitYear=COREVY,
         glu=COREGLUC, chol=CORETCHO, ldl=CORELDLC, 
         hdl=COREHDLC, tg=CORETRI) %>%
  select(subjID, visitYear, glu, chol, ldl, hdl, tg)

drawData2_whi_c1 <- read_tsv("../data/whi/phen/draws2_c1.txt", skip=10,
                             col_types=cols_only(DRAWID="c", DRAWVTYP="i", DRAWVY="i"))
drawData2_whi_c2 <- read_tsv("../data/whi/phen/draws2_c2.txt", skip=10,
                             col_types=cols_only(DRAWID="c", DRAWVTYP="i", DRAWVY="i"))
drawData2_whi <- bind_rows(drawData2_whi_c1, drawData2_whi_c2) %>%
  # filter(DRAWVTYP==1) %>%  # Only care about blood samples from first year
  filter(DRAWVTYP!=4) %>%  # Throw out non-routine (non-annual) visits
  rename(visitYear=DRAWVY) %>%
  select(DRAWID, visitYear)
labData2_whi_c1 <- read_tsv("../data/whi/phen/labs2_c1.txt", skip=10)
labData2_whi_c2 <- read_tsv("../data/whi/phen/labs2_c2.txt", skip=10)
labData2_whi <- bind_rows(labData2_whi_c1, labData2_whi_c2) %>%
  filter(SPECTYPE=="Serum",
         TESTABBR %in% c("LDLC","HDLC","TCHO","TRI","GLUC","CRP","INSU")) %>%
  inner_join(drawData2_whi, by="DRAWID") %>%
  group_by(SUBJID, visitYear, TESTABBR) %>%  # When multiple draws from the same visit...
  summarise(TESTVAL=mean(TESTVAL, na.rm=T)) %>%  # ...take the mean
  spread(key=TESTABBR, value=TESTVAL) %>%
  rename(subjID=SUBJID, chol=TCHO, ldl=LDLC, hdl=HDLC, tg=TRI, glu=GLUC, hscrp=CRP, ins=INSU)

labData_whi <- bind_rows(labData1_whi, labData2_whi) %>%
  group_by(subjID, visitYear) %>%
  summarise_all(mean, na.rm=T)  # Take mean when there are duplicate individuals from CORE and non-CORE

medsData_whi_c1 <- read_tsv("../data/whi/phen/medications_c1.txt", skip=10)
medsData_whi_c2 <- read_tsv("../data/whi/phen/medications_c2.txt", skip=10)
medsRef_whi <- read_tsv("../data/whi/phen/medication_classes.dat")
medsData_whi <- bind_rows(medsData_whi_c1, medsData_whi_c2) %>%
  # filter(F44VY==1) %>%
  inner_join(medsRef_whi, by="TCCODE") %>%
  mutate(ht_med=grepl(paste("DIURETIC", "CALCIUM BLOCKER", "ACE INHIBITOR", 
                            "ANGIOTENSIN II", "BETA BLOCKER", "BETA-BLOCKER", 
                            "ALPHA 1", "ALPHA-2", "VASODILATOR", "ALDOSTERONE",
                            collapse="|"), TCNAME),
         lipid_med=grepl("HMG COA REDUCTASE", TCNAME),
         dm_med=grepl(paste("INSULIN", "GLUCOSIDASE", "BIGUANIDE", "MEGLITINIDE", 
                      "SULFONYLUREA", "THIAZOLIDINEDIONES", collapse="|"), 
                      TCNAME)) %>%
  rename(subjID=SUBJID, visitYear=F44VY) %>%
  group_by(subjID, visitYear) %>%
  summarise(ht_med=any(ht_med),
            lipid_med=any(lipid_med),
            dm_med=any(dm_med)) %>%
  ungroup() %>%
  select(subjID, visitYear, ht_med, lipid_med, dm_med)
retroMedsData_whi <- medsData_whi %>%  
  filter(visitYear == 1) %>%
  mutate(visitYear=0)  # Assumes that meds info from visit 1 also applied at baseline
medsData_whi <- bind_rows(medsData_whi, retroMedsData_whi) %>%
  distinct()

phenoData_whi <- basicData_whi %>%
  left_join(randomization_whi, by="subjID") %>%
  left_join(behaviorData_whi, by="subjID") %>%
  left_join(Reduce(function(x,y) full_join(x,y,by=c("subjID","visitYear")),
                   list(examData_whi, labData_whi, medsData_whi)), 
            by="subjID") %>%
  mutate(subjID=as.character(subjID))

# Diet

ffq_nutrients_whi_c1 <- read_tsv("../data/whi/diet/ffq_nutrients_c1.txt", skip=10)
ffq_nutrients_whi_c2 <- read_tsv("../data/whi/diet/ffq_nutrients_c2.txt", skip=10)
ffq_nutrients_whi <- bind_rows(ffq_nutrients_whi_c1, ffq_nutrients_whi_c2) %>%
  rename(visitYear=F60VY,
         tot_cal=F60ENRGY,
         tot_fat=F60FAT,
         sfa=F60SFA,
         mufa=F60MFA,
         pufa=F60PFA,
         carb=F60CARB,
         sugars=F60TSUGR,
         sfapct=F60SFPCT,
         sucrose=F60SUCR,
         palmitic=F60SF160,
         linoleic=F60PF182,
         n3=F60OMGA3,
         n6=F60OMGA6,
         sodium=F60SODUM) %>%
  mutate(subjID=as.character(SUBJID)) %>%
  select(subjID, visitYear, tot_cal, tot_fat, sfa, mufa, pufa, carb, sugars,
         sucrose, palmitic, linoleic, n3, n6, sodium) %>%
  mutate_at(vars(-one_of("subjID","visitYear","tot_cal")), funs(./tot_cal)) %>%
  mutate(sfa2pufa=sfa/pufa,
         n62n3=n6/n3)

ffq_items_whi_c1 <- read_tsv("../data/whi/diet/ffq_items_c1.txt", skip=10)
ffq_items_whi_c2 <- read_tsv("../data/whi/diet/ffq_items_c2.txt", skip=10)
ffq_items_whi <- bind_rows(ffq_items_whi_c1, ffq_items_whi_c2) %>%
  rename(visitYear=F60VY,
         fruit=FRUITS,
         vegetables=VEGTABLS,
         fish=FISH,
         redmeat=REDMEAT,
         nuts=NUTS,
         dairy=DAIRY,
         grains=GRAINS,
         wholegrains=WHLGRNS) %>%
  mutate(subjID=as.character(SUBJID)) %>%
  select(subjID, visitYear, fruit, vegetables, fish, redmeat, 
         nuts, dairy, grains, wholegrains)

dietData_whi <- inner_join(ffq_nutrients_whi, ffq_items_whi, 
                           by=c("subjID","visitYear"))

whiMetaData <- inner_join(phenoData_whi, dietData_whi, 
                          by=c("subjID","visitYear"))
saveRDS(whiMetaData, "../int/metaData_whi.rds")
write_csv(whiMetaData, "../int/metaData_whi.csv")

# ## WHI
# whi_event_names <- c("CHD","CREVASC","STROKE")
# whi_event_days <- paste0(whi_event_names, "DY")
# outcomeData_whi_c1 <- read_tsv("../data/whi/phen/outcomes_ctos_c1.txt", skip=10) %>%
#   select(SUBJID, whi_event_names, whi_event_days, ENDFOLLOWDY)
# outcomeData_whi_c2 <- read_tsv("../data/whi/phen/outcomes_ctos_c2.txt", skip=10) %>%
#   select(SUBJID, whi_event_names, whi_event_days, ENDFOLLOWDY)
# outcomeData_whi <- bind_rows(outcomeData_whi_c1, outcomeData_whi_c2) %>%
#   filter(SUBJID %in% sampleData_whi$subjID) %>%
#   mutate(event=rowSums(.[,whi_event_names], na.rm=T)>0,
#          time=do.call(pmin, c(.[,c(whi_event_days,"ENDFOLLOWDY")], na.rm=T)),  # Event time or censoring time
#          eventType=case_when(time==CHDDY ~ "chd",
#                              time==CREVASCDY ~ "crevasc",
#                              time==STROKEDY ~ "stroke",
#                              TRUE ~ as.character(NA)),
#          incCHD=CHD==1,
#          incStroke=STROKE==1,
#          pastEvent=F) %>%
#   mutate(subjID=as.character(SUBJID)) %>%
#   select(subjID, pastEvent, event, eventType, time, incCHD, incStroke)




# ### OUTPUTS FOR DOWNSTREAM USE
# metaData <- Reduce(function(x,y) inner_join(x, y),  # Default natural join by subjID and study
#                    list(sampleData, phenoData, outcomeData))
# saveRDS(metaData, file="../int/metaData.rds")
# write_csv(metaData, "../int/metaData.csv")
# 
# sampleSheet <- inner_join(sampleData, select(metaData, subjID, sex, age, pastEvent, event), by="subjID") %>%
#   mutate(study=ifelse(grepl("lbc", study), "lbc", study))
# saveRDS(sampleSheet, file="../int/sampleSheet.rds")
