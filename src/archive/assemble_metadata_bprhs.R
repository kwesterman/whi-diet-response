library(tidyverse)
library(haven)


## Phenotypes

phenos_list <- read_sas("../data/bprhs/phen/csbl_26sep2013_list.sas7bdat")
phenos <- read_sas("../data/bprhs/phen/csbl_26sep2013.sas7bdat") %>%
  mutate(sex=ifelse(female == 1, "F", "M")) %>%
  select(studyid, sex, age, BMI, ldl, vldl, hdl, trig, sysbp, gluc, crp, mhmg) %>%
  rename(bmi=BMI,
         tg=trig,
         glu=gluc,
         hscrp=crp,
         lipid_med=mhmg)

## FFQ

ffq_bprhs <- read_sas("../data/bprhs/diet/cphhdv1b_datot_111512.sas7bdat") %>%
  rename(tot_cal=kcal_f,
         tot_fat=fat_f,
         sfa=sfa_f,
         mufa=mfa_f,
         pufa=pfa_f,
         carb=tcho_f,
         sugars=tsugar_f,
         sucrose=sucr_f,
         palmitic=s16_0_f,
         linoleic=p18_2_f,
         n3=omega3_fs,
         sodium=na_f) %>%
  mutate(n6=linoleic) %>%   # THIS IS ONLY A ROUGH APPROXIMATION FOR THE MOMENT!
  select(studyid, tot_cal, tot_fat, sfa, mufa, pufa, carb, sugars,
         sucrose, palmitic, linoleic, n3, n6, sodium) %>%
  mutate_at(vars(-one_of("studyid","tot_cal")), funs(./tot_cal)) %>%
  mutate(sfa2pufa=sfa/pufa,
         n62n3=n6/n3)


bprhsMetaData <- inner_join(phenos, ffq_bprhs, by="studyid")
write_csv(bprhsMetaData, "../int/metaData_bprhs.csv")
