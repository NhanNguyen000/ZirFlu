library(tidyverse)

ZirFlu <- list()

ZirFlu$donorInfo <- rawDat$metadata %>% 
  select(-c(probenID, Time)) %>% distinct() %>%
  select(-matches("H1N1|H3N2|B_")) %>%
  rename("patientID" = "Patient.ID")

ZirFlu$donorSamples <- rawDat$metadata %>% 
  select(Patient.ID, probenID, Time) %>% distinct() %>%
  rename("patientID" = "Patient.ID", "time" = "Time")

# antibody measure -----------------------------------------------------------------
ZirFlu$HAItiter_Feb2022 <- rawDat$HAItiter_Feb2022 %>% 
  rename("condition" = "Condition", "patientID" = "PatientID", "category" = "Category") %>%
  get.log2() %>%
  mutate(condition = tolower(condition)) %>%
  mutate(disease = ifelse(condition == "healthy", "healthy", "cirrhosis"))

ZirFlu$HAItiter_Aug2022 <- rawDat$HAItiter_Aug2022 %>%
  mutate_at(c("H1N1_d55.75", "H3N2_d55.75", "Bvictoria.Maryland_d55.75", 
              "Byamagata.Phuket_d55.75"), as.numeric) %>%
  mutate(H1N1_titerFC = ifelse(H1N1_d0 < 10, H1N1_d21.35, H1N1_d21.35/H1N1_d0),
         H3N2_titerFC = ifelse(H3N2_d0 < 10, H3N2_d21.35, H3N2_d21.35/H3N2_d0),
         Bvictoria.Maryland_titerFC = ifelse(Bvictoria.Maryland_d0 < 10, Bvictoria.Maryland_d21.35, 
                                             Bvictoria.Maryland_d21.35/Bvictoria.Maryland_d0),
         Byamagata.Phuket_titerFC = ifelse(Byamagata.Phuket_d0 < 10, Byamagata.Phuket_d21.35,
                                           Byamagata.Phuket_d21.35/Byamagata.Phuket_d0)) %>%
  rename("condition" = "Condition", "patientID" = "PatientID", "category" = "Category") %>%
  get.log2() %>%
  mutate(condition = tolower(condition)) %>%
  mutate(disease = ifelse(condition == "healthy", "healthy", "cirrhosis"))

# cytokine data -----------------------------------------------------------------
ZirFlu$cytokine_Feb2022 <- rawDat$cytokine_Feb2022 %>%
  mutate(PatientID = substring(PatientID, 1, 11)) %>%
  rename("condition" = "X...Condition", "time" = "Measure", "patientID" = "PatientID" ) %>%
  mutate(condition = tolower(condition)) %>%
  mutate(disease = ifelse(condition == "healthy", "healthy", "cirrhosis"))

# protein data -----------------------------------------------------------------
ZirFlu$protein_dat <- rawDat$protein %>%
  mutate(SampleID = stringr::str_remove(SampleID, "^0+")) %>% 
  filter(SampleID %in% ZirFlu$donorSamples$probenID) %>% 
  filter(Assay_Warning == 'PASS') %>% 
  filter(QC_Warning == 'PASS') %>% 
  filter(MissingFreq < 0.30) %>% 
  reshape2::dcast(data = ., SampleID ~ OlinkID, value.var = 'NPX') %>% 
  arrange(match(SampleID, ZirFlu$donorSamples$probenID)) %>% 
  tibble::column_to_rownames(var = 'SampleID')
all(ZirFlu$donorSamples$probenID == rownames(ZirFlu$protein_dat)) # TRUE = the order of samples are the same

ZirFlu$protein_annot <- rawDat$protein %>% 
  distinct(Assay, .keep_all = TRUE) %>% select(OlinkID, UniProt, Assay) %>%
  filter(OlinkID %in% colnames(ZirFlu$protein_dat)) %>% 
  arrange(match(OlinkID, colnames(ZirFlu$protein_dat)))
all(colnames(ZirFlu$protein_dat) == ZirFlu$protein_annot$OlinkID) # TRUE = the order of the OlinkIDs are the same

length(intersect(ZirFlu$donorSamples$probenID, colnames(rawDat$metabolite))) # 140 samples overlap, missing 18 samples

# metabolite data --------------------------------------------------------------
source("02_filterMetabolite_v2.R") # run the script to filter drug- abd food-related metabolites
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

# testing
# a <- rawDat$metabolite  %>% tibble::column_to_rownames('ionIdx') %>% 
#   select(-all_of(c("ionMz", "ionAverageInt", "ionTopFormula", "ionTopIon", "ionTopName"))) %>%
#   t() %>% as.data.frame %>% get.log2()
# a2 <- t(a); b2 <- t(dat_temp)
# all.equal(a2 %>% as_tibble %>% select('339156196'), b2 %>% as_tibble %>% select('339156157'))
# cor(a2 %>% as_tibble %>% select("339151941" ), b2 %>% as_tibble %>% select("339151988" ))

identical(rawDat$metabolite$ionIdx, as.numeric(rownames(rawDat$metabolite))) # TRUE, the ionIdx = the row index
identical(unique(ZirFlu$metabolite_annot$ionIdx), as.numeric(colnames(ZirFlu$metabolite_dat))) # TRUE, the ionIdx = the row index

setwd("/vol/projects/CIIM/Influenza/ZirrFlu/ZirFlu_NhanNguyen")
save(ZirFlu, file = "ZirFlu.RData")
