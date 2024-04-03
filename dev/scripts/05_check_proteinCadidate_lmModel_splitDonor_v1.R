library(tidyverse)

# function ----------------------
get.scatter_plot5 <- function(plotDat) {
  plotDat %>% ggplot(aes(x = abFoldchange, y = expression)) +
    geom_point(aes(color = category), size = 3, alpha = 0.8,
               position = position_jitter(w = 0.08, h = 0.1)) +
    geom_smooth(method = "lm", se = FALSE) + theme_bw()
}

# plot protein data-------------------------------
inputDat <- ZirFlu$donorInfo %>% full_join(ZirFlu$donorSamples) %>%
  select(patientID, probenID, sex, age, condition, category, time) %>%
  full_join(ZirFlu$HAItiter_Feb2022 %>% select(-condition)) %>%
  right_join(ZirFlu$protein_dat %>% rownames_to_column("probenID") %>% 
               mutate(probenID = as.numeric(probenID)))

selectedVal <- "OID20576"
abTiter <- "H1N1_titerFC"

downReg_proteins <- c("FIS1", "CCL22", "CRHBP", "TNFSF10", "MGLL", "DBNL", "TRIM21")
upReg_proteins <- c("FASLG", "TNFSF12", "CXCL8", "AGRP", "IL17D", 
                   "NCR1", "LILRB4", "ESM1", "IL15", "NTF3", "SIT1", "CEACAM21",
                   "FGF5", "MMP10", "CLEC7A", "TNFRSF13C", "CCL23", "TNFRSF13B")
check_protein <- ZirFlu$protein_annot %>% 
  filter(Assay %in% c(downReg_proteins, upReg_proteins))

# Tcell_proteins <- c("IL4R", "IL5RA", "IL10", "CCL17", "CCL22", "IL6", "IL12B")
# check_protein <- ZirFlu$protein_annot %>%
#   filter(Assay %in% Tcell_proteins)

plotList_total <- list()
for (selectedVal in check_protein$OlinkID) {
  plotList <- list()
  for(abType in c("H1N1_titerFC", "H3N2_titerFC", "Bvictoria.Maryland_titerFC", "Byamagata.Phuket_titerFC")) {
    plot_dat <- inputDat %>% filter(time == "T1") %>% 
      select(patientID, condition, disease, category, time, 
             all_of(c(abType, selectedVal))) %>%
      rename("abFoldchange" = 6, "expression" = 7)
    plotList[[abType]] <- cowplot::plot_grid(
      get.scatter_plot5(plot_dat %>% filter(disease == "healthy")),
      get.scatter_plot5(plot_dat %>% filter(disease == "cirrhosis")),
      nrow = 2
    )
  }
  plotList_total[[selectedVal]] <- plotList
}

# downReg proteins
protein <- "OID20533"
protein <- "OID20611" # maybe, need to ask
protein <- "OID20681" # not sure, need to ask
protein <- "OID20711" # good - MGLL
protein <- "OID20722" # good - FIS1
protein <- "OID20747" # not sure, need to ask, CRHBP
protein <- "OID20765" # not sure, need to ask
cowplot::plot_grid(plotList_total[[protein]]$H1N1_titerFC, 
                   plotList_total[[protein]]$H3N2_titerFC,
                   plotList_total[[protein]]$Bvictoria.Maryland_titerFC,
                   plotList_total[[protein]]$Byamagata.Phuket_titerFC,
                   nrow = 1)

# upReg_proteins
check_protein <- ZirFlu$protein_annot %>% filter(Assay %in% upReg_proteins)

protein <- "OID20429" # not sure, need to ask - but 3 strain are good
protein <- "OID20462" # maybe SIT1
protein <- "OID20480" # could be use, TNFRSF13C
protein <- "OID20481" # not sure, need to ask - but 3 strain are good, IL17D
protein <- "OID20490" # not sure
protein <- "OID20548" # could be use, CEACAM21
protein <- "OID20562" # IL15, not sure - in H3N2?
protein <- "OID20566" # could be use, not sure, NCR1
protein <- "OID20591" # need to ask, good in 3 strain, NTF3
protein <- "OID20624" # need to ask, good in 3 strain, TNFSF12
protein <- "OID20631"
protein <- "OID20636" # need to ask, good in 3 strain, CLEC7A
protein <- "OID20658" # maybe?, 	AGRP
protein <- "OID20665" # FASLG
protein <- "OID20687"
protein <- "OID20693" 	# need to ask, good in 3 strain, CCL23
protein <- "OID20702"
protein <- "OID20758" # maybe, ESM1

protein <- "OID20576" # "IL4R"
protein <- "OID20601" # "IL5RA"
protein <- "OID20431" # "IL10"
protein <- "OID20563" # "IL6"
protein <- "OID20666" # "IL12B"
protein <- "OID20745" # "CCL17", 
protein <- "OID20765" # "CCL22"
cowplot::plot_grid(plotList_total[[protein]]$H1N1_titerFC, 
                   plotList_total[[protein]]$H3N2_titerFC,
                   plotList_total[[protein]]$Bvictoria.Maryland_titerFC,
                   plotList_total[[protein]]$Byamagata.Phuket_titerFC,
                   nrow = 1)

save(plotList_total, file = "temp/plotList_total.RData")
load("temp/plotList_total.RData")

sigProteinT1 <- ZirFlu$protein_annot %>%
  mutate(sigProtein_Abpvalue = ifelse(OlinkID %in% sigProteinT1_Ab_pvalue, 1, 0),
         sigProtein_Abpadj= ifelse(OlinkID %in% sigProteinT1_Ab_padj, 1, 0),
         sigProtein_onlyAbpvalue= ifelse(OlinkID %in% sigProteinT1_onlyAb_pvalue, 1, 0),
         sigProtein_onlyAbpadj= ifelse(OlinkID %in% sigProteinT1_onlyAb_padj, 1, 0))

write.table(sigProteinT1, "temp/ZirFlu_proteinOlink_T1.txt", 
            row.names = FALSE, quote = FALSE)
write.csv(sigProteinT1, "temp/ZirFlu_proteinOlink_T1.csv", 
            row.names = FALSE, quote = FALSE)
# plot metabolite data-------------------------------
inputDat <- ZirFlu$donorInfo %>% full_join(ZirFlu$donorSamples) %>%
  select(patientID, probenID, sex, age, condition, category, time) %>%
  full_join(ZirFlu$HAItiter_Feb2022 %>% select(-condition)) %>%
  right_join(ZirFlu$metabolite_dat %>% rownames_to_column("probenID") %>% 
               mutate(probenID = as.numeric(probenID)))

metabolites <- c(2, 19, 45, 125, 190, 203, 216, 237, 240, 270, 294, 
                 302, 362, 372, 467, 668, 730, 753, 814, 877, 1166)

check_metabolite <- ZirFlu$metabolite_annot %>% 
  filter(ionIdx %in% metabolites)

plotList_total <- list()
for (selectedVal in check_metabolite$ionIdx) {
  plotList <- list()
  for(abType in c("H1N1_titerFC", "H3N2_titerFC", "Bvictoria.Maryland_titerFC", "Byamagata.Phuket_titerFC")) {
    plot_dat <- inputDat %>% filter(time == "T1") %>% 
      select(patientID, condition, disease, category, time, 
             all_of(c(abType, selectedVal))) %>%
      rename("abFoldchange" = 6, "expression" = 7)
    plotList[[abType]] <- cowplot::plot_grid(
      get.scatter_plot5(plot_dat %>% filter(disease == "healthy")),
      get.scatter_plot5(plot_dat %>% filter(disease == "cirrhosis")),
      nrow = 2
    )
  }
  plotList_total[[as.character(selectedVal)]] <- plotList
}

#save(plotList_total, file = "temp/plotList_total_metabolites.RData")
load("temp/plotList_total_metabolites.RData")

metabolite <- "2"
metabolite <- "19"
metabolite <- "45"
metabolite <- "203"
metabolite <- "216"
metabolite <- "237"
metabolite <- "294"
metabolite <- "302"
metabolite <- "362"
metabolite <- "467"
metabolite <- "730"
metabolite <- "753"
metabolite <- "814"
metabolite <- "877"
metabolite <- "125"

metabolite <- "190"
metabolite <- "240"
metabolite <- "270"
metabolite <- "668"
metabolite <- "1166"
cowplot::plot_grid(plotList_total[[metabolite]]$H1N1_titerFC, 
                   plotList_total[[metabolite]]$H3N2_titerFC,
                   plotList_total[[metabolite]]$Bvictoria.Maryland_titerFC,
                   plotList_total[[metabolite]]$Byamagata.Phuket_titerFC,
                   nrow = 1)


# CHECK Imed DATA ------
load("/vol/projects/CIIM/Influenza/ZirrFlu/iMED_NhanNguyen/iMED.RData")
# protein
inputDat_iMED <- iMED$donorInfo %>% full_join(iMED$donorSample) %>%
  select(patientID, name, gender, age, responder, time) %>% 
  full_join(iMED$HAItiter) %>% 
  right_join(iMED$protein_dat %>% rownames_to_column("name"))

selected_proteins <- c("FIS1", "CCL22", "CRHBP", "TNFSF10", "MGLL", "DBNL", "TRIM21", 
                       "FASLG", "TNFSF12", "CXCL8", "AGRP", "IL17D", "IL4R", 
                       "NCR1", "LILRB4", "ESM1", "IL15", "NTF3", "SIT1", "CEACAM21",
                       "FGF5", "MMP10", "CLEC7A", "TNFRSF13C", "CCL23", "TNFRSF13B")

check_protein <- iMED$protein_annot %>% filter(Assay %in% selected_proteins)

plotList_total <- list()
for (selectedVal in check_protein$OlinkID) {
  plotList <- list()
  for(abType in c("ab_H1N1", "ab_H3N2", "ab_B")) {
    plot_dat <- inputDat_iMED %>% filter(time == "T1") %>% 
      select(patientID, time, responder, all_of(c(abType, selectedVal))) %>%
      rename("abFoldchange" = 4, "expression" = 5, "category" = responder)
    plotList[[abType]] <- cowplot::plot_grid(
      get.scatter_plot5(plot_dat %>% filter(category != "TR")),
      get.scatter_plot5(plot_dat %>% filter(category == "TR")),
      nrow = 2
    )
  }
  plotList_total[[as.character(selectedVal)]] <- plotList
}

protein <- "OID20462" # SIT1
protein <- "OID20665" # FASLG
protein <- "OID20481" # IL17D
protein <- "OID20722" # FIS1
protein <- "OID20624" # TNFSF12
protein <- "OID20747" # CRHBP
protein <- "OID20658" # AGRP
cowplot::plot_grid(plotList_total[[protein]]$ab_H1N1, 
                   plotList_total[[protein]]$ab_H3N2,
                   plotList_total[[protein]]$ab_B,
                   nrow = 1)


# metabolites
inputDat_iMED <- iMED$donorInfo %>% full_join(iMED$donorSample) %>%
  select(patientID, name, gender, age, responder, time) %>% 
  full_join(iMED$HAItiter) %>% 
  right_join(iMED$metabolite %>% rownames_to_column("name"))


ionIdx_crossCohorts <- ZirFlu$metabolite_annot %>% filter(ionIdx %in% metabolites) %>% 
  select(ionIdx, Formula) %>% distinct() %>% rename("ionIdx_ZirFlu" = ionIdx) %>%
  left_join(iMED$metabolite_annot %>% select(ionIdx, Formula) %>% distinct()) %>%
  rename("ionIdx_iMED" = ionIdx)
check_ionIdx <- ionIdx_crossCohorts$ionIdx_iMED[!is.na(ionIdx_crossCohorts$ionIdx_iMED)]

plotList_total <- list()
for (selectedVal in check_ionIdx) {
  plotList <- list()
  for(abType in c("ab_H1N1", "ab_H3N2", "ab_B")) {
    plot_dat <- inputDat_iMED %>% filter(time == "T1") %>% 
      select(patientID, time, responder, all_of(c(abType, selectedVal))) %>%
      rename("abFoldchange" = 4, "expression" = 5, "category" = responder)
    plotList[[abType]] <- cowplot::plot_grid(
      get.scatter_plot5(plot_dat %>% filter(category != "TR")),
      get.scatter_plot5(plot_dat %>% filter(category == "TR")),
      nrow = 2
    )
  }
  plotList_total[[as.character(selectedVal)]] <- plotList
}

metabolite <- "3" # C3H6O
metabolite <- "43" # C4H6O2
metabolite <- "84" # C4H9NO2
metabolite <- "205" # C4H8N2O3
metabolite <- "290" # C8H8O3
metabolite <- "310" # C9H16O2
metabolite <- "354" #C6H12O5
metabolite <- "435" # C9H9NO3
metabolite <- "442" #C6H6N4O3
metabolite <- "515" # C6H12O7, C7H8N4O3
metabolite <- "669" # C15H24O
metabolite <- "1072" # C16H30O4
metabolite <- "1127" # C18H32O3
metabolite <- "1218" # C18H32O4
metabolite <- "1334" # C18H36O5
cowplot::plot_grid(plotList_total[[metabolite]]$ab_H1N1, 
                   plotList_total[[metabolite]]$ab_H3N2,
                   plotList_total[[metabolite]]$ab_B,
                   nrow = 1)

#save(plotList_total, file = "temp/plotList_total_metabolites_iMED.RData")
