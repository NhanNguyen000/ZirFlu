library(tidyverse)

metabolite_annot <- rawDat$metabolite_annot %>% fill(ionIdx, .direction = "down")
# HMDB endogenous metabolites -------------------
hmdb_endogenous <- read_csv("reference/20221011_HMDB_endogenousMetabolites")
annot <- metabolite_annot %>% 
  slice(which(Formula %in% hmdb_endogenous$CHEMICAL_FORMULA))
length(unique(annot$ionIdx)) # 786 metabolites

### AnnotationHub to get drugbank IDs ----------------------------------
# library(AnnotationHub)
# ah_database <- AnnotationHub()[["AH79817"]] %>% # the original ID mapping containing 9 different ID formats
#   full_join(AnnotationHub()[["AH83115"]]) %>% # data that includes common names for each compound
#   full_join(AnnotationHub()[["AH91792"]]) # current version of the mapping table that also accounts for tautomers
# save(ah_database, file = "../reference/ah_database.RData")

load("reference/ah_database.RData")
ah_info <- ah_database %>% slice(which(KEGG %in% annot$CompoundID)) %>%
  rbind(ah_database %>% slice(which(ChEBI %in% substring(annot$CompoundID, 7, 11)))) %>%
  rbind(ah_database %>% slice(which(HMDB %in% annot$CompoundID))) %>%
  distinct()

metabolite_DBids <- ah_info %>% filter(!is.na(Drugbank)) # 1345 drugbank IDs for metabolites

### DrugBank ------------------------------------------------------------------
load("reference/drugs_db.RData")
#unique(drugs_db$group)

real_drugs <- drugs_db %>% filter(group != "nutraceutical") %>% # remove drug with nutrient purpose
  filter(commercial_drug == TRUE) # only use drug have been in the market

drugs_inDat <- real_drugs %>% filter(`drugbank-id` %in% metabolite_DBids$Drugbank) # but still have nutrition function, and the drug are also natural endogenous metabolites 
drugs_inDat2 <- drugs_inDat %>% left_join(ah_info, by = c("drugbank-id" = "Drugbank")) 
drugs_inDat_list <- unique(c(drugs_inDat2$KEGG, drugs_inDat2$ChEBI, drugs_inDat2$HMDB))

annot2 <- annot %>% mutate(CompoundID = gsub("CHEBI:", "", CompoundID)) %>%
  filter(CompoundID %in% drugs_inDat_list) #68 metabolites after drug filtering steps

annot_final <- annot # use the endogenous metabolites, 786 metabolites

# remove unneeded variables ----------------------------------------------------
rm(metabolite_annot, hmdb_endogenous, metabolite_DBids, ah_database, ah_info, 
   real_drugs, drugs_db, drugs_inDat, annot, annot2 )
