# WGCNA per time point ==================================================================
library("WGCNA")
options(stringsAsFactors = FALSE)

## Protein data - for T1 --------------------------------------------
proteinDat <- ZirFlu$proteinImputeDat_time2$T1
gsg <- goodSamplesGenes(proteinDat, verbose = 3) # check the data to be sure
gsg$allOK

sampleTree <- get.WGCNA_sampleTree(proteinDat)
get.WGCNA_powerTables(proteinDat)
WGCNA_network <- get.WGCNAnet(proteinDat, selected_power = 4, netType = "unsigned")

table(WGCNA_network$colors) # the modules and theirs size, the lable 0 is for genes outside of all modules
WGCNA_modules <- get.module_elements(WGCNA_network)
MEs <- orderMEs(moduleEigengenes(proteinDat, 
                                 labels2colors(WGCNA_network$colors))$eigengenes)
datTraits <- ZirFlu$log2_abTiters_Feb2022 %>%
  select(matches("titerFC|d0|d21.35"), PatientID) %>% 
  left_join(ZirFlu$metadata %>% select(probenID, Patient.ID), 
            by = c("PatientID" = "Patient.ID")) %>%
  filter(probenID %in% rownames(MEs)) %>% column_to_rownames("probenID") %>%
  select(-PatientID) %>%
  rename_with(~gsub("Byamagata.Phuket", "Byamagata", .x, fixed = TRUE)) %>%
  rename_with(~gsub("Bvictoria.Maryland", "Bvictoria", .x, fixed = TRUE))

WGCNA_moduleTrait <- get.moduleTrait_relation(MEs, datTraits)

# between control and cirrhosis
a<- MEs %>% rownames_to_column("probenID") %>% 
  mutate(probenID = as.numeric(probenID)) %>%
  left_join(ZirFlu$metadata)
a %>% ggplot(aes(x = condition, y = MEturquoise)) + 
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1))

library(ggpubr)
ggboxplot(a, x = "condition", y = "MEturquoise", 
          color = "condition", palette = "jco", add = "jitter") + 
  stat_compare_means(method = "t.test",
                     comparisons = list(c("healthy", "decompensated cirrhosis"), 
                                        c("decompensated cirrhosis", "compensated cirrhosis"),
                                        c("healthy", "compensated cirrhosis"))) +
  stat_compare_means(method = "anova", label.y = 0.5) 

# find top proteins
WGCNA_modules_annot <- list()
for (name in names(WGCNA_modules)) {
  WGCNA_modules_annot[[name]] <- get.proteinAnnot(ZirFlu$protein_annot, WGCNA_modules[[name]])
}

chooseTopHubInEachModule(proteinDat, WGCNA_network$colors)
connectivityTable <- intramodularConnectivity.fromExpr(proteinDat, 
                                                       colors =  WGCNA_network$colors)
rownames(connectivityTable) <- names(proteinDat)
intraMod_connect <- get.intraModule_connectivity(proteinDat,
                                                 WGCNA_network, WGCNA_modules_annot)
intraMod_connect$turquoise %>% arrange(desc(kWithin)) %>% 
  filter(kWithin >2.6) %>% select(Assay)
# check protein TNFRSF14 (OID20783)
a2 <- proteinDat %>% rownames_to_column("probenID") %>% 
  mutate(probenID = as.numeric(probenID)) %>%
  left_join(ZirFlu$metadata)
ggboxplot(a2, x = "condition", y = "OID20783", 
          color = "condition", palette = "jco", add = "jitter") + 
  stat_compare_means(method = "t.test",
                     comparisons = list(c("healthy", "decompensated cirrhosis"), 
                                        c("decompensated cirrhosis", "compensated cirrhosis"),
                                        c("healthy", "compensated cirrhosis"))) +
  stat_compare_means(method = "anova", label.y = 0.5) 

# correlation
a3 <- a2 %>% 
  left_join(ZirFlu$log2_abTiters_Feb2022, by = c("Patient.ID" = "PatientID"))

a3 %>% 
  ggplot(aes(x = OID20783, y = H1N1_d21.35)) + 
  geom_point(aes(shape = Condition, color = Category)) +
  geom_smooth(method = "lm", se = FALSE) +
  xlab("TNFRSF14 (OID20783) level at T1 (baseline)")

library(ggpmisc)
a3 %>% 
  ggplot(aes(x = OID20783, y = H1N1_d21.35)) + 
  geom_point(aes(shape = Condition, color = Category)) +
  stat_poly_line(se = FALSE) + 
  stat_poly_eq(aes(label = paste(after_stat(eq.label),
                                 after_stat(rr.label), sep = "*\", \"*")))+
  xlab("TNFRSF14 (OID20783) level at T1 (baseline)")