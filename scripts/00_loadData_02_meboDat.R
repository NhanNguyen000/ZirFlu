rm(list = ls())

library(readxl)
library(tidyverse)

# functions used in the analysis -----------------------------------------------
get.log2 <- function(dat) {
  library(tidyverse)
  library(rstatix)
  outcome <- dat %>% mutate(across(where(is.numeric), ~log2(.x))) %>%
    mutate(across(where(is.numeric), ~ifelse(is.infinite(.x), 0, .x)))
  return(outcome)
}

# load the data list object ----------------------------------------------------
load("data/ZirFlu.RData")

# load metabolite data -------------------------------------------------------
metabolite <- read_excel(
  '/vol/projects/BIIM/Influenza/ZirrFlu/metabolic/raw_data/spreadsheets/DATA.xlsx',
  sheet = 'ions')

metabolite_annot <- read_excel(
  '/vol/projects/BIIM/Influenza/ZirrFlu/metabolic/raw_data/spreadsheets/DATA.xlsx',
  sheet = 'annotation') %>% fill(ionIdx, .direction = "down")

# filter drug- and food-related metabolites ------------------------------------

## HMDB endogenous metabolites ----------------------------------------------
hmdb_endogenous <- read_csv("reference/20221011_HMDB_endogenousMetabolites")
endoMebo_annot <- metabolite_annot %>% slice(which(Formula %in% hmdb_endogenous$CHEMICAL_FORMULA))
#length(unique(endoMebo_annot$ionIdx)) # 786 metabolites

## check whether metabolites relate to drug--------------------------------------

# ### get drugbank IDs from AnnotationHub
# # library(AnnotationHub)
# # ah_database <- AnnotationHub()[["AH79817"]] %>% # the original ID mapping containing 9 different ID formats
# #   full_join(AnnotationHub()[["AH83115"]]) %>% # data that includes common names for each compound
# #   full_join(AnnotationHub()[["AH91792"]]) # current version of the mapping table that also accounts for tautomers
# # save(ah_database, file = "reference/ah_database.RData")
# 
# load("reference/ah_database.RData") # get drugbank IDs per metabolite from AnnotationHub
# ah_info <- ah_database %>% slice(which(KEGG %in% endoMebo_annot$CompoundID)) %>%
#   rbind(ah_database %>% slice(which(ChEBI %in% substring(endoMebo_annot$CompoundID, 7, 11)))) %>%
#   rbind(ah_database %>% slice(which(HMDB %in% endoMebo_annot$CompoundID))) %>%
#   distinct()
# 
# metabolite_DBids <- ah_info %>% filter(!is.na(Drugbank)) # 1345 drugbank IDs for metabolites
# 
# ### get drug information from DrugBank 
# load("reference/drugs_db.RData") # extract from DrugBank database xml file in 2022
# 
# real_drugs <- drugs_db %>% filter(group != "nutraceutical") %>% # remove drug with nutrient purpose
#   filter(commercial_drug == TRUE) # only use drug have been in the market
# 
# ### drug-related metabolites in both AnnotationHub and DrugBank 
# drugs_inDat <- real_drugs %>% filter(`drugbank-id` %in% metabolite_DBids$Drugbank) # but still have nutrition function, and the drug are also natural endogenous metabolites 
# drugs_inDat2 <- drugs_inDat %>% left_join(ah_info, by = c("drugbank-id" = "Drugbank")) 
# drugs_inDat_list <- unique(c(drugs_inDat2$KEGG, drugs_inDat2$ChEBI, drugs_inDat2$HMDB))

# check if endogenous metabolites are drug-related metabolites in both AnnotationHub and DrugBank
# annot_withoutDrug <- endoMebo_annot %>% mutate(CompoundID = gsub("CHEBI:", "", CompoundID)) %>%
#   filter(CompoundID %in% drugs_inDat_list) # 216 metabolites after drug filtering steps
# 
# rm(metabolite_annot, hmdb_endogenous, metabolite_DBids, ah_database, ah_info,
#    real_drugs, drugs_db, drugs_inDat, drugs_inDat2, annot_withoutDrug )

## metabolite annotation used----------------------------------------------------
ZirFlu$metabolite_annot <- endoMebo_annot

# sample IDs correction ----------------------------------------------------
dat_temp <- metabolite  %>% tibble::column_to_rownames('ionIdx') %>% 
  select(-all_of(c("ionMz", "ionAverageInt", "ionTopFormula", "ionTopIon", "ionTopName"))) %>%
  t() %>% as.data.frame %>% get.log2()

old_probenID <- c(339151941, 339156196, 339151948, 339156278, 339151926, 339152850,
                  339156227, 339156287, 339156214, 339156220, 339156213, 339156314,
                  339156322, 339152826, 339152815, 339152794, 339156212, 339151960)
correct_probenID <- c(339151988, 339156157, 339151947, 339156279, 339151924, 339152857,
                      339156180, 339156286, 339156199, 339156204, 339156219, 339156288,
                      339156321, 339152802, 339152806, 339152807, 339156181, 339152000)

change_probenID <- data.frame("old_probenID" = old_probenID, "correct_probenID" = correct_probenID)
for (id in 1:nrow(change_probenID)) {
  rownames(dat_temp)[which(rownames(dat_temp) == change_probenID$old_probenID[id])] <- change_probenID$correct_probenID[id]
}

ZirFlu$metabolite_dat <-  dat_temp[which(rownames(dat_temp) %in% ZirFlu$donorSamples$probenID),] %>%
  select(unique(ZirFlu$metabolite_annot$ionIdx))

identical(metabolite$ionIdx, as.numeric(rownames(metabolite))) # TRUE, the ionIdx = the row index
identical(unique(ZirFlu$metabolite_annot$ionIdx), as.numeric(colnames(ZirFlu$metabolite_dat))) # TRUE, the ionIdx = the row index

rm(endoMebo_annot, dat_temp, old_probenID, correct_probenID, change_probenID, id)

## save data -------------------------------------------------------------------
save(ZirFlu, file = "data/ZirFlu.RData")
