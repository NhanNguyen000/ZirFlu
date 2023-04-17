library("KEGGREST")
library(tidyverse)

### Annotate using formula ----------------------------------

# KEGG database: using formula to annotate metabolite (able to annotate half of the metabolite formula)
inputDat <- ZirFlu$metabolite_annot %>% 
  select(ionIdx, Formula) %>% distinct() # 799 formulas (786 Ids)

kegg_annot <- c()
for (formula in inputDat$Formula) {
  annot_temp <- keggFind("compound", formula, "formula")
  exact_formula <- annot_temp[which(annot_temp == formula)]
  kegg_annot <- c(kegg_annot, exact_formula)
}

kegg_outcomes <- kegg_annot %>% enframe() %>% distinct() %>%
  rename("cpdId" = "name", "Formula" = "value") %>% 
  left_join(inputDat) # 524 formulas (516 ionIdx)

keggID_allMetabolites <- kegg_outcomes %>% mutate(cpdId = str_replace(cpdId, "cpd:", ""))
save(keggID_allMetabolites, 
     file = "temp/20221222_keggID_allMetabolites.RData")

### get kegg pathway database ----------------------------------
cpdpathway <- keggLink("cpd", "pathway") %>% enframe %>%
  rename("mapId" = "name", "cpdId" = "value")
listpathway <- keggList("pathway") %>% enframe %>%
  rename("mapId" = "name", "pathway" = "value")
cpdpathway_combine <- cpdpathway %>% full_join(listpathway)

keggID_pathway <- kegg_outcomes %>% 
  inner_join(cpdpathway_combine) # consist of 375 formulas with 368 Ids / total 799 formulas with 786 Ids

keggID_libraryFormula <- keggID_pathway %>%
  select(pathway, Formula) %>% group_by(pathway) %>%
  summarise(across(everything(), str_c, collapse = "; ")) # 288 pathway in total

keggID_libraryIonId <- keggID_pathway %>%
  select(pathway, ionIdx) %>% group_by(pathway) %>%
  summarise(across(everything(), str_c, collapse = "; ")) # 288 pathway in total

save(keggID_libraryFormula, keggID_libraryIonId, 
     file = "temp/20221019_keggIDusingFormula_library.RData")

# check kegg Brite -----------------------------
keggIDs_temp <- kegg_outcomes %>% as.data.frame() %>%
  separate(col = "cpdId", into = c("cpd", "keggId"), sep = ":")
kegg_IDs <- keggIDs_temp$keggId

kegg_db <- map_dfr(keggGet(kegg_IDs[1:(1+9)]), ~as_tibble(t(.)))
for (i in seq(11, length(kegg_IDs), 10)) {
  kegg_db <- kegg_db %>%
    full_join(map_dfr(keggGet(kegg_IDs[i:(i+9)]), ~as_tibble(t(.))))
}

kegg_brite <- get.keggBrite(kegg_db)
kegg_brite_taxonomy <- get.keggTaxonomy(kegg_brite)

unique(unlist(kegg_brite)[grep("BR:br", unlist(kegg_brite))]) # 20 different groups
kegg_brite_v2 <- get.kegg_briteGroups(kegg_brite)

