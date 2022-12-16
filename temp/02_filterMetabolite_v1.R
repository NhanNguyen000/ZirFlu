library(tidyverse)
# library(dplyr)
# library(openxlsx)
# library(magrittr)
# library(ggplot2)
# library(reshape2)
# library(ggsci)

# need to check Maritjn code?

# drug metabolite ----------------------------
# drugIDs_info <- read.csv2("reference/20220907_drugBankIDs_info_NN.csv") # manual check in DrugBank database, need to update
# drugs_butAre_humanMetabolites <- c("Testosterone", "Dehydroepiandrosterone",
#                                    "9-cis Retinoic acid", "17-Methyltestosterone",
#                                    "Acetate", "Phenylacetic acid")
# 
# real_drugs <- drugIDs_info %>% filter(HaveBrandNames == "Yes") %>% # select drug that are on the market
#   filter(InGroupNutraceutical != "Yes") %>% # select drug that is not for nutrient purpose
#   filter(Note == "") %>% slice(-c(2, 26, 29, 35)) %>% # remove duplicated drugbank IDs 
#   slice(-which(Name %in% drugs_butAre_humanMetabolites)) # keep drug are also human metabolites
# 
# drug_notInDrugBank <- c("aprobarbital", "acetazolamide", "paracetamol", 
#                        "amoxicilin", "propericiazine")
# metabolies_noDrug <- ZirFlu$metabolite_annot %>% 
#   slice(-which(Formula %in% real_drugs$X...Top.annotation.formula)) %>%
#   slice(-grep(paste(drug_notInDrugBank, collapse = "|"), tolower(CompoundName)))
# 
# length(unique(ZirFlu$metabolite_annot$ionIdx))
# length(unique(metabolies_noDrug$ionIdx))

# HMDB endogenous metabolites -------------------
hmdb_endogenous <- read_csv("reference/20221011_HMDB_endogenousMetabolites")
annot <- ZirFlu$metabolite_annot %>% 
  slice(which(Formula %in% hmdb_endogenous$CHEMICAL_FORMULA))
length(unique(annot$ionIdx))

# annot2 <- annot %>% 
#   slice(-which(Formula %in% real_drugs$X...Top.annotation.formula)) %>%
#   slice(-grep(paste(drug_notInDrugBank, collapse = "|"), tolower(CompoundName)))
# length(unique(annot2$ionIdx))

### AnnotationHub to get drugbank IDs ----------------------------------
library(AnnotationHub)
ah_database <- AnnotationHub()[["AH79817"]] %>% # the original ID mapping containing 9 different ID formats
  full_join(AnnotationHub()[["AH83115"]]) %>% # data that includes common names for each compound
  full_join(AnnotationHub()[["AH91792"]]) # current version of the mapping table that also accounts for tautomers
save(ah_database, file = "temp/ah_database.RData")

load("temp/ah_database.RData")
ah_info <- ah_database %>% slice(which(KEGG %in% annot$CompoundID)) %>%
  rbind(ah_database %>% slice(which(ChEBI %in% substring(annot$CompoundID, 7, 11)))) %>%
  rbind(ah_database %>% slice(which(HMDB %in% annot$CompoundID))) %>%
  distinct()

metabolite_DBids <- ah_info %>% filter(!is.na(Drugbank)) # 441 metabolites

### DrugBank ----------------------------------
load("reference/drugs_db.RData")
#unique(drugs_db$group)

real_drugs <- drugs_db %>% filter(group != "nutraceutical") %>% # remove drug with nutrient purpose
  filter(commercial_drug == TRUE) # only use drug have been in the market
drugs_inDat <- real_drugs %>% filter("drugbank-id" %in% ah_info$Drugbank)
dim(drugs_inDat) # nrow = 0, no drugs in the metabolites data

annot_final <- annot
