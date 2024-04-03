# get sample at baseline
sample_T1 <- ZirFlu$donorSamples %>% filter(time == "T1")

# proteins
selected_proteins <- ZirFlu$protein_annot %>% 
  filter(Assay %in% c("TNFSF10", "CXCL8", "IL6", "IL17C", "IL17D", "TNFSF12"))

protein_dat <- ZirFlu$protein_dat %>% select(selected_proteins$OlinkID) %>% 
  rownames_to_column("probenID") %>% filter(probenID %in% sample_T1$probenID)

# metabolites
metabolite_dat <- ZirFlu$metabolite_dat %>% 
  select(unique(c(metaSig_onlyAb_2019, metaSig_onlyAb_2020))) %>%
  rownames_to_column("probenID") %>% filter(probenID %in% sample_T1$probenID)

metabolite_dat <- ZirFlu$metabolite_dat %>% 
  select(selected_metabolites) %>% # metabolite 75% consistent across samples & season
  rownames_to_column("probenID") %>% filter(probenID %in% sample_T1$probenID)

# check the association
library(rstatix)

outcome <- c()
for (protein in selected_proteins$OlinkID) {
  outcome_temp <- c()
  
  for (metabolite in names(metabolite_dat)[-1]) {
    dat_temp <- protein_dat %>% select(probenID, any_of(protein)) %>%
      full_join(metabolite_dat %>% select(c("probenID", any_of(metabolite)))) %>%
      column_to_rownames("probenID")
    
    cor_res <- cor_test(dat_temp)
    outcome_temp <- rbind(outcome_temp, cor_res)
  }
  outcome[[protein]] <- outcome_temp
}

outcome_pval <- outcome %>% lapply(function(x) x %>% filter(p < 0.05))

g3 <- ZirFlu$metabolite_annot %>% 
  filter(ionIdx %in% outcome_pval$OID20611$var2) %>% 
  select(ionIdx, Formula) %>% distinct()


g4 <- ZirFlu$metabolite_annot %>% 
  filter(ionIdx %in% outcome_pval$OID20631$var2) %>% 
  select(ionIdx, Formula) %>% distinct()

outcome_all <- outcome_pval %>% 
  reduce(full_join) %>% mutate(var2 = as.numeric(var2)) %>% 
  left_join(ZirFlu$protein_annot %>% select(OlinkID, Assay),
            by = c("var1" = "OlinkID")) %>% 
  left_join(ZirFlu$metabolite_annot %>% select(ionIdx, Formula) %>% distinct(), 
            by = c("var2" = "ionIdx")) %>%
  relocate(c(Assay, Formula)) %>% select(-c(var1, var2))

outcome_all  %>% 
  ggplot(aes(x = Assay, y = Formula)) + 
  geom_tile(aes(fill = cor)) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") + theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# only select metabolites show consistent trend (>75% across strains and years)
metabolites_allPro <- outcome_all %>% 
  add_count(Formula) %>% filter(n == 6) %>% select(-n)

metabolites_allPro  %>% 
  ggplot(aes(x = Assay, y = Formula)) + 
  geom_tile(aes(fill = cor)) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") + theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
 

load("temp/20221222_keggID_allMetabolites.RData")
g4 <- keggID_allMetabolites %>% filter(ionIdx %in% as.numeric(g2$var2))

