library(OlinkAnalyze)
library(readxl)
library(tidyverse)

# 23.02.2023 - after talking with Saumya and Marthij about Olink protein data from plasma
# realized that protein T2 (for for visit 1 at day 3-7 after vaccine) and T3 (for visit 2 at day 21-28 after vaccine) 
# differ from abTiter T2 (for visit 2 at day 28 after vaccine) and T3 (for visit 3 at day >50/60 after vaccine)

setwd("/vol/projects/CIIM/Influenza/ZirrFlu/ZirFlu_NhanNguyen")

proteinPlate <- read_excel("../proteomic/raw_data/ZirrFlu plates final.xlsx",
                           sheet = "Phenotypes")
proteinSample <- proteinPlate %>% select(-1, -c(HLA_HEP, Date, Tissue)) %>%
  rename("sex" = "Sex", "age" = "Age",  "condition" = "Condition", 
         "time" = "Time", "season" = "Season") %>%
  filter(season != "Bridge") 

# proteinSample_plan <- read_excel("../proteomic/raw_data/Auszug-ZirFlu 1 und 2 CentraXX23112021_YangLi.xlsx",
#                                       sheet = "Tabelle1") %>%
#   mutate(probenID = stringr::str_remove(ProbenID, "^0+")) %>% 
#   filter(probenID %in% ZirFlu$donorSamples$probenID) %>%
#   select(ProbandenID, GeburtsJahr, Geschlecht, DatumProbe, Therapiephase, Zeitpunkt, probenID) %>% 
#   distinct()
# proteinSample_test <- ZirFlu$donorSamples %>% full_join(proteinSample_plan)
# identical(proteinSample_test$patientID, proteinSample_test$ProbandenID) # TRUE

ZirFlu$donorSamples <- proteinSample %>%
  select(patientID, probenID, season, time) %>%
  mutate(probenID = stringr::str_remove(probenID, "^0+")) %>%
  mutate(time = ifelse(time == "Baseline", "T1", 
                       ifelse(time == "T1", "T2", "T3")))

ZirFlu$donorInfo2 <- proteinSample %>%
  select(patientID, sex, age, condition) %>% distinct()

ZirFlu$protein_dat <- rawDat$protein %>%
  slice(which(SampleID %in% proteinSample$probenID)) %>% 
  filter(Assay_Warning == 'PASS') %>% 
  filter(QC_Warning == 'PASS') %>% 
  filter(MissingFreq < 0.30) %>% 
  mutate(SampleID = stringr::str_remove(SampleID, "^0+")) %>% 
  reshape2::dcast(data = ., SampleID ~ OlinkID, value.var = 'NPX') %>% 
  arrange(match(SampleID, donorSamples$probenID)) %>% 
  tibble::column_to_rownames(var = 'SampleID')
all(ZirFlu$donorSamples$probenID == rownames(ZirFlu$protein_dat)) # TRUE = the order of samples are the same
