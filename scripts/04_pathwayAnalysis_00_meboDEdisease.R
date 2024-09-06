rm(list = ls())

library(tidyverse)

# Aim of this code:
# this code will generate files that contains metabolite/compound IDs, 
# the list of metabolite IDs can be used to run pathway analysis in https://www.metaboanalyst.ca/ 
# the list of metabolite IDs consists of HMDB, CHEBI, and KEGG ids, but it can be used directly in the website
# You can select the input type HMDB or KEGG ids, the result will be the same
# Note: https://www.metaboanalyst.ca/ does not have option to add the background metabolites to avoid potential bias in the pathway analysis
# however, we still use this website because so far it has the most comprehensive metabolite databases to do pathway analysis 

# load the data------------------------------------------------------
load("processedData/ZirFlu.RData")

# DE metabolte associated to disease -----------------------------------------------------------------
# outcome of the lm() model between metabolite and health conditions (cirrhosis vs. healthy), 
load("processedData/lmRes_disease.RData")
load("processedData/DEmebo_disease.RData")

meboDEdisease_2019 <- venn_disease_2019 %>% flatten() %>% unlist()
length(unique(meboDEdisease_2019))

meboDEdisease_annot <- ZirFlu$metabolite_annot %>% 
  filter(Formula %in% unique(meboDEdisease_2019))

# make the metabolite/compound ID file
write.table(meboDEdisease_annot$CompoundID, 
            file = "output/meboDEdiseaseAnnot_season2019.txt", sep = "\t",
            row.names = FALSE, col.names = FALSE, quote = FALSE)

# DE metabolte associated to ab titer -----------------------------------------------------------------
# outcome of the lm() model between antibody titer visit2 and metabolites
load("processedData/lmRes_meboAbVisit2.RData")
load("processedData/DEmebo_abVisit2.RData")

meboDEabTiter_2019 <- venn_meboAbVisit2_2019 %>% flatten() %>% unlist()
length(unique(meboDEabTiter_2019))

meboDEabTiter_annot <- ZirFlu$metabolite_annot %>% 
  filter(Formula %in% unique(meboDEabTiter_2019))

# make the metabolite/compound ID file
write.table(meboDEabTiter_annot$CompoundID, 
            file = "output/meboDEabTiterAnnot_season2019.txt", sep = "\t",
            row.names = FALSE, col.names = FALSE, quote = FALSE)

# DE metabolte associated to ab titer, consistent across strains and seasons -------------------------------------
load("processedData/DEmebo_abVisit2_consistAcrossStrainSeason.RData")

meboDEabTiter_consistent <- ZirFlu$metabolite_annot %>% 
  filter(Formula %in% unique(selected_lmStatistic$Formula))

# make the metabolite/compound ID file
write.table(meboDEabTiter_consistent$CompoundID, 
            file = "output/meboDEabTiterAnnot_season2019_consistent.txt", sep = "\t",
            row.names = FALSE, col.names = FALSE, quote = FALSE)
