library(tidyverse)

### WHI ###

whi_items_c1 <- read_tsv("../data/whi/diet/ffq_items_c1.txt", skip=10)
whi_items_c2 <- read_tsv("../data/whi/diet/ffq_items_c2.txt", skip=10)
whi_items <- bind_rows(whi_items_c1, whi_items_c2) %>%
  filter(F60VTYP==1,  # Only use FFQ from first visit
         STATUS==3) %>%  # Status of 1 or 2 indicates suspect FFQ quality (unreasonable total energy)
  select(SUBJID, FRUITS:DAIRY)  # Variables "FRUITS" through "DAIRY" are computed food group levels

whi_nutrients_c1 <- read_tsv("../data/whi/diet/ffq_nutrients_c1.txt", skip=10)
whi_nutrients_c2 <- read_tsv("../data/whi/diet/ffq_nutrients_c2.txt", skip=10)
whi_nutrients <- bind_rows(whi_nutrients_c1, whi_nutrients_c2) %>%
  filter(F60VTYP==1,
         STATUS==3)


### FHS ###

fhs_ffq_c1 <- read_tsv("../data/fhs/diet/ffq_c1.txt", skip=10)
fhs_ffq_c2 <- read_tsv("../data/fhs/diet/ffq_c2.txt", skip=10)
fhs_ffq <- bind_rows(fhs_ffq_c1, fhs_ffq_c2)
