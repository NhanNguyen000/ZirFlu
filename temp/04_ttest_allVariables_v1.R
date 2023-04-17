library(rstatix)
meta_ttest_outcome <- list()
abTiter_ttest <- ZirFlu$metadata %>% 
  full_join(ZirFlu$abTiters_Feb2022, by = c("Patient.ID" = "PatientID")) %>%
  select(-probenID) %>% distinct() %>%
  mutate(group = ifelse(condition == "healthy", "healthy", "cirrhosis")) %>%
  get.log2()

meta_ttest_outcome[["abTiter"]] <- get.ttest(dat = abTiter_ttest, 
                                             selected_var = names(abTiter_ttest)[14:21])
a <- ZirFlu$metadata %>% 
  full_join(ZirFlu$log2_abTiters_Feb2022, by = c("Patient.ID" = "PatientID")) %>%
  select(-probenID) %>% distinct() %>%
  mutate(group = ifelse(condition == "healthy", "healthy", "cirrhosis"))
a2 <- list()
a2[["abTiter"]] <- get.ttest(dat = a, selected_var = names(a)[14:21])

#my_data_ttest %>% group_by(group) %>% get_summary_stats(H1N1_d0, type = "mean_sd")

cytokine_ttest <- ZirFlu$metadata %>% 
  full_join(ZirFlu$cytokine, by = c("Patient.ID" = "PatientID")) %>%
  mutate(across(where(is.numeric), ~na_if(., 0.01))) %>% replace(is.na(.), 0) %>%
  select(-probenID) %>% distinct() %>%
  mutate(group = ifelse(condition == "healthy", "healthy", "cirrhosis")) %>%
  get.log2()
meta_ttest_outcome[["cytokine_d0"]] <- get.ttest(cytokine_ttest %>% filter(Measure == "d0"), 
                                                 selected_var = names(cytokine_ttest)[16:48])
meta_ttest_outcome[["cytokine_d21-35"]] <- get.ttest(cytokine_ttest %>% filter(Measure == "d21-35"), 
                                                     selected_var = names(cytokine_ttest)[16:48])

a2[["cytokine_d0"]] <- get.ttest(ZirFlu$log2_cytokine %>% filter(Measure == "d0") %>% 
                                   rename(Patient.ID = PatientID, group = Condition), 
                                 selected_var = names(ZirFlu$log2_cytokine)[5:37])
abTiter_pairedTtest <- ZirFlu$metadata %>% 
  full_join(ZirFlu$abTiters_Feb2022, by = c("Patient.ID" = "PatientID")) %>%
  select(-probenID) %>% distinct() %>%
  mutate(condition = ifelse(condition == "healthy", "healthy", "cirrhosis")) %>%
  get.log2()

outcome_temp <- as.data.frame(matrix(nrow = 0, ncol = 9))
names(abTiter_pairedTtest_outcome ) <- c("condition", ".y.", "group1", "group2", "n1", "n2","statistic", "df", "p")
for (i in c("H1N1", "H3N2", "Bvictoria.Maryland", "Byamagata.Phuket")) {
  dat_temp <- abTiter_pairedTtest %>%
    gather(key = "group", value = {{i}}, paste0(i, "_d0"), paste0(i, "_d21.35")) %>%
    select(Patient.ID, group, {{i}}, condition) %>% rename(value = {{i}})
  res_temp <- dat_temp %>% group_by(condition) %>% t_test(value ~ group, paired = TRUE)
  res_temp[,2] <- i
  outcome_temp  <- rbind(outcome_temp, res_temp)
}

meta_ttest_outcome[["abTiter_pairedTtest"]] <- outcome_temp

# cytokine 
meta_ttest_outcome[["cytokine_pairedTtest"]] <- get.pairedTtest(cytokine_ttest %>% mutate(group = Measure) %>%
                                                                  mutate(condition = ifelse(condition == "healthy", "healthy", "cirrhosis")),
                                                                names(cytokine_ttest)[16:48])
# For protein data - no complete case: 62 at T1, 37 at T3, 59 at T4
protein_ttest <- merge(ZirFlu$metadata2, ZirFlu$protein_dat, by.x = "probenID", by.y = "row.names")
# protein ttest time T1 vs T3
protein_T1_T3 <- get.completeCase_dat(protein_ttest, time_point = c("T1", "T3")) %>%
  rename(group = Time)

protein_ttest_T1_T3 <- list()
protein_ttest_T1_T3[["total"]] <- get.ttest(dat = protein_T1_T3,
                                            selected_var = names(protein_T1_T3)[30:337])
protein_ttest_T1_T3[["healthy"]] <- get.ttest(protein_T1_T3 %>% filter(condition == "healthy"), 
                                              selected_var = names(protein_T1_T3)[30:337])
protein_T1_T3_cirrhosis_dat <- protein_T1_T3 %>% filter(condition != "healthy")
# some proteins do not have enough observation to perform ttest
rm_proteins <- c("OID20533", "OID20614", "OID20617", "OID20644",
                 "OID20646", "OID20656", "OID20666", "OID20668",
                 "OID20709", "OID20778")
#View(protein_T1_T3_cirrhosis_dat %>% 
#       select(c(c("Patient.ID", "condition", "group"), rm_proteins)))
possible_proteins <- setdiff(names(protein_T1_T3)[30:337], rm_proteins)
protein_ttest_T1_T3[["cirrhosis"]] <- get.ttest(protein_T1_T3_cirrhosis_dat,
                                                selected_var = possible_proteins)
# protein ttest T1 vs. T4
protein_T1_T4 <- get.completeCase_dat(protein_ttest, time_point = c("T1", "T4")) %>%
  rename(group = Time)
protein_ttest_T1_T4 <- list()
protein_ttest_T1_T4[["total"]] <- get.ttest(dat = protein_T1_T4, 
                                            selected_var = names(protein_T1_T4)[30:337])
protein_ttest_T1_T4[["healthy"]] <- get.ttest(protein_T1_T4 %>% filter(condition == "healthy"), 
                                              selected_var = names(protein_T1_T4)[30:337])
protein_ttest_T1_T4[["cirrhosis"]] <- get.ttest(protein_T1_T4 %>% filter(condition != "healthy"), 
                                                selected_var = names(protein_T1_T4)[30:337])

# paired ttest - protein T1 vs. T3

protein_pairedTtest <- list()
rm_proteins <- c("OID20533", "OID20614", "OID20617", "OID20644",
                 "OID20646", "OID20656", "OID20666", "OID20668",
                 "OID20709", "OID20752", "OID20778")
#View(protein_T1_T3 %>% 
#       select(c(c("Patient.ID", "condition", "group"), rm_proteins)))
possible_proteins <- setdiff(names(protein_T1_T3)[30:337], rm_proteins)

protein_pairedTtest[["T1_T3"]] <- get.pairedTtest(protein_T1_T3 %>%
                                                    mutate(condition = ifelse(condition == "healthy", "healthy", "cirrhosis")),
                                                  possible_proteins)
protein_pairedTtest[["T1_T4"]] <- get.pairedTtest(protein_T1_T4 %>% 
                                                    mutate(condition = ifelse(condition == "healthy", "healthy", "cirrhosis")),
                                                  names(protein_T1_T4)[30:337])
# check proteins with p < 0.05
meta_ttest_outcome$abTiter %>% filter(p<0.05)
meta_ttest_outcome$cytokine_d0 %>% filter(p<0.05) 
meta_ttest_outcome$`cytokine_d21-35` %>% filter(p<0.05)
meta_ttest_outcome$abTiter_pairedTtest %>% filter(p<0.05)
meta_ttest_outcome$cytokine_pairedTtest %>% filter(p<0.05) # 3 cytokines

View(protein_ttest_T1_T3$total %>% filter(p<0.05)) # 27 proteins
View(protein_ttest_T1_T3$healthy %>% filter(p<0.05)) # 31 proteins
View(protein_ttest_T1_T3$cirrhosis %>% filter(p<0.05)) # 1 protein - OID20507 - not see at T4

View(protein_ttest_T1_T4$total %>% filter(p<0.05)) # 15 proteins
View(protein_ttest_T1_T4$healthy %>% filter(p<0.05)) # 29 proteins
View(protein_ttest_T1_T4$cirrhosis %>% filter(p<0.05)) # no protein

View(protein_pairedTtest$T1_T3 %>% filter(p<0.05)) # 37 proteins, only healthy peoples
View(protein_pairedTtest$T1_T4 %>% filter(p<0.05)) # 32 proteins, only healthy peoples

# For metabolite data - no complete case: 44 at T1, 37 at T3, 59 at T4
metabolite_ttest <- merge(ZirFlu$metadata2, ZirFlu$metabolite_dat, by.x = "probenID", by.y = "row.names")
View(metabolite_ttest[which(metabolite_ttest$Time == "T3"),]) # mostly, healthy patients
View(metabolite_ttest[which(metabolite_ttest$Time == "T4"),])
