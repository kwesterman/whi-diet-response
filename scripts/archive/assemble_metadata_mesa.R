library(tidyverse)


## Phenotypes

phenos <- read_tsv("../data/mesa/phen/phs000209.v10.pht001116.v7.p2.c1.MESA_Exam1Main.GRU.txt",
                   skip=10) %>% 
  rename(subjID=sidno,
         age=age1c,
         sex=gender1,
         race=race1c,
         bmi=bmi1c,
         smk_now=cig1c,
         smk_py=pkyrs1c,
         diabetes_recode=dm031t,
         hypertension=htn1c,
         ldl=ldl1,
         hdl=hdl1,
         chol=chol1,
         tg=trig1,
         glu=glucos1c,
         sbp=sbp1c,
         crp=crp1,
         statin=sttn1c) %>%
  mutate(sex=ifelse(sex == 0, "F", "M")) %>%
  select(subjID, sex, age, race, bmi, smk_now, smk_py,
         diabetes_recode, hypertension, statin,
         ldl,  hdl, chol, tg, glu, sbp, crp)
  
  
## FFQ

ffq <- read_tsv("../data/mesa/diet/phs000209.v10.pht002107.v3.p2.c1.MESA_Exam1DietNutrients.GRU.txt",
                skip=10) %>%
  rename(subjID=sidno,
         sfa=tsfan1c,
         mufa=tmufan1c,
         pufa=tpufan1c,
         sfa_pct=pclsfn1c,
         mufa_pct=pclmfn1c,
         pufa_pct=pclpfn1c) %>%
  select(subjID, sfa, mufa, pufa, sfa_pct, mufa_pct, pufa_pct)


mesa_meta_data <- inner_join(phenos, ffq, by="subjID")
write_csv(mesa_meta_data, "../int/metaData_mesa.csv")
