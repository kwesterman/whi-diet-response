library(tidyverse)
library(kinship2)

ped_df <- read_tsv("../data/fhs/fhs_pedigree.txt", skip=10)
ped <- with(ped_df, pedigree(shareid, dadid=fshare, momid=mshare, 
                             sex=sex))
kmat <- kinship(ped)
kmat <- kmat[as.character(ped$id), as.character(ped$id)]

write_csv(ped, "../../int/fhs_pedigree.csv")
saveRDS(kmat, "../../int/fhs_kinship_matrix.rds")
