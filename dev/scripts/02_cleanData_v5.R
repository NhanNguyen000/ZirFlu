library(tidyverse)

setwd("/vol/projects/CIIM/Influenza/ZirrFlu/ZirFlu_NhanNguyen")
ZirFlu <- list()
## HAI antibody titers ---------------------------------------------------------
ZirFlu$HAItiter_2019 <- rawDat$HAItiter_2019 %>% 
  mutate(condition = tolower(condition),
         vaccine_response = tolower(vaccine_response),
         category = ifelse(vaccine_response == "non", "NR",
                           ifelse(vaccine_response %in% c("single", "double"), "Other", "TR"))) %>%
  relocate(condition, .after = last_col()) %>%
  rename_with(~sub("d0", "T1", .x)) %>%
  rename_with(~sub("d21-35", "T2", .x)) %>% 
  rename_with(~sub("d55-75", "T3", .x)) %>%
  mutate_at(c(2:17), as.numeric) %>% get.log2()

ZirFlu$HAItiter_2020 <- rawDat$HAItiter_2020 %>%
  mutate(condition = tolower(condition),
         vaccine_response = tolower(vaccine_response),
         category = ifelse(vaccine_response == "non", "NR",
                           ifelse(vaccine_response %in% c("single", "double"), "Other", "TR"))) %>%
  relocate(condition, .after = last_col()) %>%
  rename_with(~sub("d0", "T1", .x)) %>%
  rename_with(~sub("d21-35", "T2", .x)) %>% 
  rename_with(~sub("d56-104", "T3", .x)) %>%
  mutate_at(c(2:17), as.numeric) %>% get.log2() %>%
  filter(patientID != "Z-01-99-069")
# The participant Z-01-99-069 in season 2020 has azathioprine in his medication 
# and should be excluded from the analysis (info from the doctor).

ZirFlu$HAItiter <- ZirFlu$HAItiter_2019 %>% mutate(season = "2019") %>%
  full_join(ZirFlu$HAItiter_2020 %>% mutate(season = "2020")) %>%
  relocate("season", .after = "patientID") %>%
  mutate(condition = factor(condition, 
                            levels = c("decompensated cirrhosis", 
                                       "compensated cirrhosis", "healthy"))) %>%
  mutate(disease = ifelse(condition == "healthy", "healthy", "cirrhosis")) %>%
  mutate(disease = factor(disease,
                          levels = c("cirrhosis", "healthy")))

## metadata --------------------------------------------------------------------
ZirFlu$donorInfo <- ZirFlu$HAItiter_2019 %>% 
  select(patientID, condition) %>% mutate(season = "2019") %>%
  full_join(ZirFlu$HAItiter_2020 %>% 
              select(patientID, condition) %>% mutate(season = "2020")) %>%
  mutate(condition = factor(condition, 
                            levels = c("decompensated cirrhosis", 
                                       "compensated cirrhosis", "healthy"))) %>%
  mutate(disease = ifelse(condition == "healthy", "healthy", "cirrhosis")) %>%
  mutate(disease = factor(disease,
                          levels = c("cirrhosis", "healthy")))

# ZirFlu$donorInfo_2019 <- rawDat$metadata_2019 %>% 
#   select(-c(probenID, Time)) %>% distinct() %>%
#   select(-matches("H1N1|H3N2|B_")) %>%
#   rename("patientID" = "Patient.ID")
# 
# ZirFlu$donorSamples_2019 <- rawDat$metadata_2019 %>% 
#   select(Patient.ID, probenID, Time) %>% distinct() %>%
#   rename("patientID" = "Patient.ID", "time" = "Time")

# update donor info based on metadata of protein samples -----------------------
proteinPlate <- read_excel("../proteomic/raw_data/ZirrFlu plates final.xlsx",
                           sheet = "Phenotypes")
proteinSample <- proteinPlate %>% select(-1, -c(HLA_HEP, Date, Tissue)) %>%
  rename("sex" = "Sex", "age" = "Age",  "condition" = "Condition", 
         "time" = "Time", "season" = "Season") %>%
  filter(season != "Bridge") 

ZirFlu$donorInfo <- proteinSample %>% # update donor info based on metadata of protein samples
  select(patientID, sex, age, condition, season) %>% distinct() %>%
  mutate(disease = tolower(condition)) %>% select(-condition) %>%
  right_join(ZirFlu$donorInfo)

ZirFlu$donorSamples <- proteinSample %>%
  select(patientID, probenID, season, time) %>%
  mutate(probenID = stringr::str_remove(probenID, "^0+")) %>%
  mutate(time = ifelse(time == "Baseline", "T1", 
                       ifelse(time == "T1", "T2", "T3")))

# protein data -----------------------------------------------------------------
ZirFlu$protein_dat <- rawDat$protein %>%
  slice(which(SampleID %in% proteinSample$probenID)) %>% 
  filter(Assay_Warning == 'PASS') %>% 
  filter(QC_Warning == 'PASS') %>% 
  filter(MissingFreq < 0.30) %>% 
  mutate(SampleID = stringr::str_remove(SampleID, "^0+")) %>% 
  reshape2::dcast(data = ., SampleID ~ OlinkID, value.var = 'NPX') %>% 
  arrange(match(SampleID, ZirFlu$donorSamples$probenID)) %>% 
  tibble::column_to_rownames(var = 'SampleID')
all(ZirFlu$donorSamples$probenID == rownames(ZirFlu$protein_dat)) # TRUE = the order of samples are the same

ZirFlu$protein_annot <- rawDat$protein %>% 
  distinct(Assay, .keep_all = TRUE) %>% select(OlinkID, UniProt, Assay) %>%
  filter(OlinkID %in% colnames(ZirFlu$protein_dat)) %>% 
  arrange(match(OlinkID, colnames(ZirFlu$protein_dat)))
all(colnames(ZirFlu$protein_dat) == ZirFlu$protein_annot$OlinkID) # TRUE = the order of the OlinkIDs are the same

rm(proteinPlate, proteinSample)
# metabolite ----------------------------------------------------------------
setwd("scripts"); source("02_filterMetabolite_v2.R") # run the script to filter drug- abd food-related metabolites
ZirFlu$metabolite_annot <- annot_final 

dat_temp <- rawDat$metabolite  %>% tibble::column_to_rownames('ionIdx') %>% 
  select(-all_of(c("ionMz", "ionAverageInt", "ionTopFormula", "ionTopIon", "ionTopName"))) %>%
  t() %>% as.data.frame %>% get.log2()

old_probenID <- c(339151941, 339156196, 339151948, 339156278, 339151926, 339152850,
                  339156227, 339156287, 339156214, 339156220, 339156213, 339156314,
                  339156322, 339152826, 339152815, 339152794, 339156212, 339151960)
correct_probenID <- c(339151988, 339156157, 339151947, 339156279, 339151924, 339152857,
                      339156180, 339156286, 339156199, 339156204, 339156219, 339156288,
                      339156321, 339152802, 339152806, 339152807, 339156181, 339152000)
change_probenID <- data.frame("old_probenID" = old_probenID, "correct_probenID" = correct_probenID)
for (i in 1:nrow(change_probenID)) {
  rownames(dat_temp)[which(rownames(dat_temp) == change_probenID$old_probenID[i])] <- change_probenID$correct_probenID[i]
}

ZirFlu$metabolite_dat <-  dat_temp[which(rownames(dat_temp) %in% ZirFlu$donorSamples$probenID),] %>%
  select(unique(ZirFlu$metabolite_annot$ionIdx))

identical(rawDat$metabolite$ionIdx, as.numeric(rownames(rawDat$metabolite))) # TRUE, the ionIdx = the row index
identical(unique(ZirFlu$metabolite_annot$ionIdx), as.numeric(colnames(ZirFlu$metabolite_dat))) # TRUE, the ionIdx = the row index

rm(annot_final, dat_temp, old_probenID, correct_probenID, change_probenID, i)
## save data -------------------------------------------------------------------
setwd("/vol/projects/CIIM/Influenza/ZirrFlu/ZirFlu_NhanNguyen")
save(ZirFlu, file = "ZirFlu.RData")
