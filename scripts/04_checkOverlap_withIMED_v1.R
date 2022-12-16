library(tidyverse)

# function ------------------------------------
get.sig_var <- function(inputDat) {
  sig_vars <- list()
  for (i in names(inputDat)) {
    sig_vars[[i]] <- inputDat[[i]] %>% filter(p.value < 0.05)
  }
  return(sig_vars)
}
# check protein ------------------------------------
load("/vol/projects/CIIM/Influenza/iMED/proteomic/antibodyAssocAnalysis/cohort2_allProtsAssoc_EachStrainEachTimePoint.RData")
# Saumya: the order of columns is H1N1_T1,H3N2_T1, B_T1, H1N1_T3,H3N2_T3,B_T3,H1N1_T4,H3N2_T4,B_T4
dim(alltime_allStrains_assoc)
iMED_protein <- list()
for (i in seq(1, 54, 6)) {
  iMED_protein[[as.character(i)]] <- alltime_allStrains_assoc[,c(i:(i+5))]
}
names(iMED_protein) <- c("H1N1_T1", "H3N2_T1", "B_T1", 
                         "H1N1_T3", "H3N2_T3", "B_T3", 
                         "H1N1_T4", "H3N2_T4", "B_T4")

iMED_sigProteins <- get.sig_var(iMED_protein)
iMED_sigProteinT1 <- c(iMED_sigProteins$H1N1_T1$Protein, 
                       iMED_sigProteins$H3N2_T1$Protein, 
                       iMED_sigProteins$B_T1$Protein)

a <- get.proteinAnnot(ZirFlu$protein_annot, 
                      sigProteinT1_Ab_pvalue) 
intersect(a$Assay, iMED_sigProteinT1)

a <- get.proteinAnnot(ZirFlu$protein_annot, 
                      sigProteinT1_Ab_padj) 
intersect(a$Assay, iMED_sigProteinT1)

a <- get.proteinAnnot(ZirFlu$protein_annot, 
                      sigProteinT1_onlyAb_pvalue) 
intersect(a$Assay, iMED_sigProteinT1)

a <- get.proteinAnnot(ZirFlu$protein_annot, 
                      sigProteinT1_onlyAb_padj) 
intersect(a$Assay, iMED_sigProteinT1)
# check metabolites ------------------------------------
iMED_metabolite <- read.csv("/vol/projects/CIIM/Influenza/iMED/metabolic/analysis/brendan/output/assoc_metabolites_lfc.csv")
load("/vol/projects/CIIM/Influenza/ZirrFlu/iMED_NhanNguyen/iMED.RData")
View(iMED$seleted_metabolite_annot)
iMED_metabolite_T1sig <- list()
for (type in c("ab_H1N1", "ab_H3N2", "ab_B")) {
  temp_dat <- iMED_metabolite %>% 
    filter(time == "T1") %>% filter(strain == type) %>%
    filter(p.value < 0.05) %>% 
    slice(which(metabolite %in% iMED$seleted_metabolite_annot$ionIdx))
  iMED_metabolite_T1sig[[type]] <- iMED$metabolite_annot %>% 
    select(ionIdx, Formula) %>% distinct() %>%
    slice(which(ionIdx %in% temp_dat$metabolite))
}

iMED_metabolite_T1 <- unique(c(iMED_metabolite_T1sig$ab_H1N1$Formula, 
                               iMED_metabolite_T1sig$ab_H3N2$Formula, 
                               iMED_metabolite_T1sig$ab_B$Formula))
# formular annotation
get.metaboliteAnnot(ZirFlu$metabolite_annot, metaboliteT1_disease_pvalue) %>%
  select(ionIdx, Formula) %>% distinct() %>% 
  filter(Formula %in% iMED_metabolite_T1)

get.metaboliteAnnot(ZirFlu$metabolite_annot, metaboliteT1_Ab_pvalue) %>%
  select(ionIdx, Formula) %>% distinct() %>% 
  filter(Formula %in% iMED_metabolite_T1)

get.metaboliteAnnot(ZirFlu$metabolite_annot, metaboliteT1_disease_padj) %>%
  select(ionIdx, Formula) %>% distinct() %>% 
  filter(Formula %in% iMED_metabolite_T1)

get.metaboliteAnnot(ZirFlu$metabolite_annot, metaboliteT1_Ab_padj) %>%
  select(ionIdx, Formula) %>% distinct() %>% 
  filter(Formula %in% iMED_metabolite_T1)

get.metaboliteAnnot(ZirFlu$metabolite_annot, metaboliteT1_onlyAb_pvalue) %>%
  select(ionIdx, Formula) %>% distinct() %>% 
  filter(Formula %in% iMED_metabolite_T1)

get.metaboliteAnnot(ZirFlu$metabolite_annot, metaboliteT1_onlyAb_padj) %>%
  select(ionIdx, Formula) %>% distinct() %>% 
  filter(Formula %in% iMED_metabolite_T1)
