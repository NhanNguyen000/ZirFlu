# WGCNA =======================================
library("WGCNA")
options(stringsAsFactors = FALSE)

load("20220908_impute_proteinDat_zirFlu_NN.RData")
library("WGCNA")
options(stringsAsFactors = FALSE)
imputed_protein <- impute_dat_norm.predict

# check the data to be sure
gsg <- goodSamplesGenes(imputed_protein, verbose = 3)
gsg$allOK

sampleTree <- get.WGCNA_sampleTree(imputed_protein) # sample human-000025 is an outliner
imputed_protein2 <- rm.samples(sampleTree, imputed_protein, cutoff = 40)
sampleTree2 <- get.WGCNA_sampleTree(imputed_protein2)
get.WGCNA_powerTables(imputed_protein2)
WGCNA_network <- get.WGCNAnet(imputed_protein2, selected_power = 4, netType = "signed")

table(WGCNA_network$colors) # the modules and theirs size, the lable 0 is for genes outside of all modules
MEs2 <- WGCNA_network$MEs
#geneTree <- WGCNA_network$dendrograms[[1]]
WGCNA_modules <- get.module_elements(WGCNA_network)

WGCNA_modules_annot <- list()
for (name in names(WGCNA_modules)) {
  WGCNA_modules_annot[[name]] <- get.proteinAnnot(ZirFlu$protein_annot, WGCNA_modules[[name]])
}

chooseTopHubInEachModule(imputed_protein2, WGCNA_network$colors)
connectivityTable <- intramodularConnectivity.fromExpr(imputed_protein2, 
                                                       colors =  WGCNA_network$colors)
rownames(connectivityTable) <- names(imputed_protein2)

connectivity_turquoise <- connectivityTable %>% rownames_to_column("OlinkID") %>% 
  filter(OlinkID %in% WGCNA_modules_annot$turquoise$OlinkID) %>%
  full_join(WGCNA_modules_annot$turquoise)

connectivity_turquoise %>% arrange(desc(kWithin)) %>% 
  filter(kWithin >2.6) %>% select(Assay) 

connectivity_blue <- connectivityTable %>% rownames_to_column("OlinkID") %>% 
  filter(OlinkID %in% WGCNA_modules_annot$blue$OlinkID) %>%
  full_join(WGCNA_modules_annot$blue)

connectivity_blue %>% arrange(desc(kWithin)) %>% 
  filter(kWithin > 9.5) %>% select(Assay) 
