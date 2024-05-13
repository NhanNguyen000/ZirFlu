rm(list = ls())

library(readxl)
library(tidyverse)

# load the data list object------------------------------------------------------
load("data/ZirFlu.RData")

# donor information (metadata) --------------------------------------------------

## metadata for HAI titer data --------------------------------------------------
donorInfo_fromHAItiter <- ZirFlu$HAItiter %>% 
  select(patientID, condition, season, condition, disease)

## update donor info based on metadata of protein samples -----------------------
proteinPlate <- read_excel("/vol/projects/BIIM/Influenza/ZirrFlu/proteomic/raw_data/ZirrFlu plates final.xlsx",
                           sheet = "Phenotypes")

proteinSample <- proteinPlate %>% select(-1, -c(HLA_HEP, Date, Tissue)) %>%
  rename("sex" = "Sex", "age" = "Age",  "condition" = "Condition", 
         "time" = "Time", "season" = "Season") %>%
  filter(season != "Bridge") 

ZirFlu$donorInfo <- proteinSample %>% # update donor info based on metadata of protein samples
  select(patientID, sex, age, condition, season) %>% distinct() %>%
  mutate(disease = condition) %>% select(-condition) %>%
  right_join(donorInfo_fromHAItiter)

# save data -------------------------------------------------------------------
# base on the metafile for sample randomisation before send out for protein and metabolite mesurement
ZirFlu$donorSamples <- proteinSample %>%
  select(patientID, probenID, season, time) %>%
  mutate(probenID = stringr::str_remove(probenID, "^0+")) %>%
  mutate(time = ifelse(time == "Baseline", "Baseline", # Day 0 is baseline
                       ifelse(time == "T1", "Visit1", # Day 3-4 is Visit 1
                              ifelse(time == "T2", "Visit2", time)))) # day 21-28 is Visit 2

save(ZirFlu, file = "data/ZirFlu.RData")
