library(tidyverse)

ZirFlu <- list()

ZirFlu$metadata <- rawDat$metadata %>%
  select(Patient.ID, season, sex, age, HEP_HLA_B27, 
         sample_baseline_date, medication, Disease, 
         child_pugh_score, cirrhosis.ethiology, cirrhosis.therapy, 
         condition, probenID)

ZirFlu$metadata2 <- rawDat$metadata 

ZirFlu$donorsInfo <- ZirFlu$metadata %>% select(-probenID) %>% distinct() # need to change to rawDat

ZirFlu$abTiters_Feb2022 <- rawDat$HAItiter_Feb2022
ZirFlu$log2_abTiters_Feb2022 <- rawDat$HAItiter_Feb2022 %>% 
  get.log2() %>%
  mutate(Condition = ifelse(Condition == "Healthy", "Healthy", "Cirrhosis"))

ZirFlu$abTiters_Aug2022 <- rawDat$HAItiter_Aug2022 %>%
  mutate(H1N1_titerFC = ifelse(H1N1_d0 < 10, H1N1_d21.35, H1N1_d21.35/H1N1_d0),
         H3N2_titerFC = ifelse(H3N2_d0 < 10, H3N2_d21.35, H3N2_d21.35/H3N2_d0),
         Bvictoria.Maryland_titerFC = ifelse(Bvictoria.Maryland_d0 < 10, Bvictoria.Maryland_d21.35, 
                                             Bvictoria.Maryland_d21.35/Bvictoria.Maryland_d0),
         Byamagata.Phuket_titerFC = ifelse(Byamagata.Phuket_d0 < 10, Byamagata.Phuket_d21.35,
                                           Byamagata.Phuket_d21.35/Byamagata.Phuket_d0))

ZirFlu$log2_abTiters_Aug2022 <- ZirFlu$abTiters_Aug2022 %>%
  get.log2() %>% mutate(donor = Condition) %>%
  mutate(Condition = ifelse(Condition == "Healthy", "Healthy", "Cirrhosis"))

ZirFlu$cytokine <- rawDat$cytokine_Aug2022 %>%
  mutate(PatientID = substring(PatientID, 1, 11))

ZirFlu$protein_dat <- rawDat$protein %>%
  mutate(SampleID = stringr::str_remove(SampleID, "^0+")) %>% 
  filter(SampleID %in% ZirFlu$metadata$probenID) %>% 
  filter(Assay_Warning == 'PASS') %>% 
  filter(QC_Warning == 'PASS') %>% 
  filter(MissingFreq < 0.30) %>% 
  reshape2::dcast(data = ., SampleID ~ OlinkID, value.var = 'NPX') %>% 
  arrange(match(SampleID, ZirFlu$metadata$probenID)) %>% 
  tibble::column_to_rownames(var = 'SampleID')
#all(ZirFlu$metadata$probenID == rownames(ZirFlu$protein_dat)) # TRUE, if the order of samples are the same

ZirFlu$protein_annot <- rawDat$protein %>% 
  distinct(Assay, .keep_all = T) %>% select(OlinkID, UniProt, Assay) %>%
  filter(OlinkID %in% colnames(ZirFlu$protein_dat)) %>% 
  arrange(match(OlinkID, colnames(ZirFlu$protein_dat)))
#all(colnames(ZirFlu$protein_dat) == ZirFlu$protein_annot$OlinkID) # TRUE, if the order of the OlinkIDs are the same

length(intersect(ZirFlu$metadata$probenID, colnames(rawDat$metabolite))) # 140 samples overlap, missing 18 samples

ZirFlu$metabolite_dat <- rawDat$metabolite  %>% 
  select(intersect(ZirFlu$donorSamples$probenID, colnames(rawDat$metabolite))) %>%
  t() %>% as.data.frame %>% get.log2()

ZirFlu$metabolite_annot <- rawDat$metabolite %>% 
  select(ionIdx, ionMz, ionAverageInt, ionTopFormula, ionTopIon, ionTopName)
identical(rawDat$metabolite$ionIdx, as.numeric(rownames(rawDat$metabolite))) # TRUE, the ionIdx = the row index

save(ZirFlu, file = "ZirFlu.RData")
