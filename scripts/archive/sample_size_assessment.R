genotyped <- read_tsv("sample_info.txt", skip=15) %>%
  rename(SUBJID=SubjectID) %>%
  distinct(SUBJID, .keep_all=T)

ffq_c1 <- read_tsv("diet/ffq_items_c1.txt", skip=10)
ffq_c2 <- read_tsv("diet/ffq_items_c2.txt", skip=10)
ffq_multiple <- bind_rows(ffq_c1, ffq_c2) %>%
  group_by(SUBJID) %>%
  filter(n()>2) %>%
  ungroup() %>%
  distinct(SUBJID, .keep_all=T)

randomization_c1 <- read_tsv("randomization_c1.txt", skip=10)
randomization_c2 <- read_tsv("randomization_c2.txt", skip=10)
randomization_DM <- bind_rows(randomization_c1, randomization_c2) %>%
  filter(DMFLAG==1,
         DMARM==1) %>%
  distinct(SUBJID, .keep_all=T)

length(unique(genotyped$SUBJID))
length(unique(ffq_multiple$SUBJID))
length(unique(randomization_DM$SUBJID))
geno_ffq_dm <- Reduce(intersect, list(genotyped$SUBJID, 
                                      ffq_multiple$SUBJID, 
                                      randomization_DM$SUBJID))

methMeta <- readRDS("../../../meth_cvd/int/metaData.rds")
length(intersect(methMeta$subjID, geno_ffq_dm))
