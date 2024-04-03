library(OlinkAnalyze)
library(readxl)
library(tidyverse)

setwd("/vol/projects/CIIM/Influenza/ZirrFlu/ZirFlu_NhanNguyen")

proteinPlate <- read_excel("../proteomic/raw_data/ZirrFlu plates final.xlsx",
                sheet = "Phenotypes")
proteinSample <- proteinPlate %>% select(-1, -c(HLA_HEP, Date, Tissue)) %>%
  rename("sex" = "Sex", "age" = "Age",  "condition" = "Condition", 
         "time" = "Time", "season" = "Season") %>%
  filter(season != "Bridge") 

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
