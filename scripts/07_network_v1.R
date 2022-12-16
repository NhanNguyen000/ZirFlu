# https://mr.schochastics.net/material/netVizR/
#  https://briatte.github.io/ggnet/
# https://bookdown.org/markhoff/social_network_analysis/culture-and-networks.html

# try the network tool im Cytoscape
library("KEGGREST")
library(tidyverse)

inputDat <- get.metaboliteAnnot(ZirFlu$metabolite_annot, metaboliteT1_disease_pvalue) %>%
  select(ionIdx, Formula) %>% distinct()

kegg_annot <- c()
for (formula in inputDat$Formula) {
  annot_temp <- keggFind("compound", formula, "formula")
  exact_formula <- annot_temp[which(annot_temp == formula)]
  kegg_annot <- c(kegg_annot, exact_formula)
}

kegg_outcomes <- kegg_annot %>% enframe() %>% distinct() %>%
  rename("keggId" = "name", "Formula" = "value") %>% 
  mutate(keggId = substring(keggId, 5, 10)) %>%
  left_join(inputDat) # 800 keggID (266 ionIdx)

write.csv(kegg_outcomes, "temp/T1_sigMetabolites.csv", 
          quote = FALSE, row.names = FALSE)

# network option 2: 
View(highlight_metabolites$Formula) # 66 metabolite formulas

kegg_annot <- c()
for (formula in highlight_metabolites$Formula) {
  annot_temp <- keggFind("compound", formula, "formula")
  exact_formula <- annot_temp[which(annot_temp == formula)]
  kegg_annot <- c(kegg_annot, exact_formula)
}

kegg_outcomes <- kegg_annot %>% enframe() %>% distinct() %>%
  rename("keggId" = "name", "Formula" = "value") %>% 
  mutate(keggId = substring(keggId, 5, 10)) %>%
  left_join(inputDat) # 228 keggID

write.csv(kegg_outcomes, "temp/T1_66highlitedMetabolites.csv", 
          quote = FALSE, row.names = FALSE)
