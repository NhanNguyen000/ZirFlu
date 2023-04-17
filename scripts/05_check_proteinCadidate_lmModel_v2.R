inputDat <- ZirFlu$donorInfo %>% full_join(ZirFlu$HAItiter) %>%
  left_join(ZirFlu$donorSamples %>% filter(time == "T1")) %>%
  left_join(ZirFlu$protein_dat %>% rownames_to_column("probenID"))

inputDat <- ZirFlu$donorInfo %>% full_join(ZirFlu$HAItiter) %>%
  left_join(ZirFlu$donorSamples %>% filter(time == "T1")) %>%
  left_join(ZirFlu$protein_dat %>% rownames_to_column("probenID")) %>%
  filter(season == "2019")
## T1 protein plot -------------------
a <- ZirFlu$protein_annot %>% filter(Assay %in% highlight_proteins$Assay)

selectedPro <- "OID20722"
abTiter <- "H1N1_abFC"

plot_dat <- inputDat %>% filter(time == "T1") %>% 
  select(patientID, condition, disease, category, {{abTiter}},
         OID20481, OID20624, OID20765) # IL17D, TNFSF12, CCL22
cowplot::plot_grid(
  ggplot(plot_dat, aes(x = H1N1_abFC, y = OID20481)) +
    geom_point(aes(color = category, shape = condition), 
               size = 3, alpha = 0.8,
               position = position_jitter(w = 0.08, h = 0.1)) +
    geom_smooth(method = "lm", se = FALSE) + theme_bw(),
  ggplot(plot_dat, aes(x = H1N1_abFC, y = OID20624)) +
    geom_point(aes(color = category, shape = condition), 
               size = 3, alpha = 0.8,
               position = position_jitter(w = 0.08, h = 0.1)) +
    geom_smooth(method = "lm", se = FALSE) + theme_bw(),
  ggplot(plot_dat, aes(x = H1N1_abFC, y = OID20765)) +
    geom_point(aes(color = category, shape = condition), 
               size = 3, alpha = 0.8,
               position = position_jitter(w = 0.08, h = 0.1)) +
    geom_smooth(method = "lm", se = FALSE) + theme_bw(),
  nrow = 3
)

plot_dat$category <- factor(plot_dat$category, levels = c("NR", "Other", "TR"))
plot_dat$disease <- factor(plot_dat$disease, 
                           levels = c("healthy", "cirrhosis"))
plot_dat$condition <- factor(plot_dat$condition, 
                             levels = c("healthy", "compensated cirrhosis", "decompensated cirrhosis"))
cowplot::plot_grid(
  ggboxplot(plot_dat, x = "category", y = "OID20481", palette = "jco", add = "jitter"),
  ggboxplot(plot_dat, x = "disease", y = "OID20481", palette = "jco", add = "jitter"),
  ggboxplot(plot_dat, x = "condition", y = "OID20481", palette = "jco", add = "jitter") + 
    rotate_x_text(20),
  ggboxplot(plot_dat, x = "category", y = "OID20624", palette = "jco", add = "jitter"),
  ggboxplot(plot_dat, x = "disease", y = "OID20624", palette = "jco", add = "jitter"),
  ggboxplot(plot_dat, x = "condition", y = "OID20624", palette = "jco", add = "jitter") + 
    rotate_x_text(20),
  ggboxplot(plot_dat, x = "category", y = "OID20765", palette = "jco", add = "jitter"),
  ggboxplot(plot_dat, x = "disease", y = "OID20765", palette = "jco", add = "jitter"),
  ggboxplot(plot_dat, x = "condition", y = "OID20765", palette = "jco", add = "jitter") + 
    rotate_x_text(20),
  nrow = 3, byrow = FALSE
)

library(tidyverse)

# function ----------------------
get.scatter_plot5 <- function(plotDat) {
  plotDat %>% ggplot(aes(x = abFoldchange, y = expression)) +
    geom_point(aes(color = category), size = 3, alpha = 0.8,
               position = position_jitter(w = 0.08, h = 0.1)) +
    geom_smooth(method = "lm", se = FALSE) + theme_bw()
}

get.scatter_plot6 <- function(plotDat) {
  plotDat %>% ggplot(aes(x = abFoldchange, y = expression)) +
    geom_point(aes(color = category), size = 3, alpha = 0.8,
               position = position_jitter(w = 0.08, h = 0.1)) +
    geom_smooth(method = "lm", se = FALSE) + theme_bw() +
    facet_wrap(vars(disease))
}

check_protein <- get.proteinAnnot(ZirFlu$protein_annot, unlist(proSig))
check_protein <- ZirFlu$protein_annot %>% 
  filter(Assay %in% highlight_proteins$Assay)
check_protein <- ZirFlu$protein_annot %>% 
  filter(OlinkID %in% selected_proteins)

plotList_total <- list()
plotList_total2 <- list()
plotList_total3 <- list()
plotList_total4 <- list()
for (selectedVal in check_protein$OlinkID) {
  plotList <- list()
  plotList2 <- list()
  plotList3 <- list()
  plotList4 <- list()
  for(abType in c("H1N1_abFC", "H3N2_abFC", "Bvictoria_abFC", "Byamagata_abFC")) {
    plot_dat <- inputDat %>% filter(time == "T1") %>% 
      select(patientID, condition, disease, category, time, 
             all_of(c(abType, selectedVal))) %>%
      rename("abFoldchange" = 6, "expression" = 7)
    plotList[[abType]] <- cowplot::plot_grid(
      get.scatter_plot5(plot_dat %>% filter(disease == "healthy")),
      get.scatter_plot5(plot_dat %>% filter(disease == "cirrhosis")),
      nrow = 2
    )
    plotList2[[abType]] <- get.scatter_plot6(plot_dat)
    plotList2[[abType]] <- get.scatter_plot6(plot_dat)
    plotList3[[abType]] <- plot_dat %>% 
      ggscatter(x = "abFoldchange", y = "expression", 
                color = "disease", add = "reg.line") +
      stat_cor(aes(color = disease), method = "pearson")
    plotList4[[abType]] <- plot_dat %>% 
      ggscatter(x = "abFoldchange", y = "expression", add = "reg.line") +
      stat_cor(method = "pearson")
  }
  plotList_total[[selectedVal]] <- plotList
  plotList_total2[[selectedVal]] <- plotList2
  plotList_total3[[selectedVal]] <- plotList3
  plotList_total4[[selectedVal]] <- plotList4
}

protein <- "OID20462" # SIT1
protein <- "OID20480" # TNFRSF13C
protein <- "OID20490" # FGF5
protein <- "OID20494" # JCHAIN
protein <- "OID20504" # ITM2A
protein <- "OID20533" # TRIM21
protein <- "OID20547" # CLEC4C
protein <- "OID20562" # IL15
protein <- "OID20574" # OSM
protein <- "OID20596" # KRT19
protein <- "OID20611" # TNSFS10
protein <- "OID20636" # CLEC7A
protein <- "OID20641" # GMPR
protein <- "OID20642" # TPSAB1
protein <- "OID20650" # VEGFA
protein <- "OID20655" # CCL13
protein <- "OID20666" # IL12B
#selectedVal <- "OID20666"
protein <- "OID20672" # MMP1
protein <- "OID20681" # DBNL
protein <- "OID20687" # MMP10
protein <- "OID20693" # CCL23
protein <- "OID20709" # CCN2
protein <- "OID20711" # MGLL
protein <- "OID20714" # SHMT1
protein <- "OID20722" # FIS1
protein <- "OID20727" # GAL
protein <- "OID20743" # EPCAM
protein <- "OID20747" # CRHBP
protein <- "OID20753" # MEPE
cowplot::plot_grid(plotList_total3[[protein]]$H1N1_abFC, 
                   plotList_total3[[protein]]$H3N2_abFC,
                   plotList_total3[[protein]]$Bvictoria_abFC,
                   plotList_total3[[protein]]$Byamagata_abFC,
                   nrow = 2)
cowplot::plot_grid(plotList_total4[[protein]]$H1N1_abFC, 
                   plotList_total4[[protein]]$H3N2_abFC,
                   plotList_total4[[protein]]$Bvictoria_abFC,
                   plotList_total4[[protein]]$Byamagata_abFC,
                   nrow = 2)
cowplot::plot_grid(plotList_total[[protein]]$H1N1_abFC, 
                   plotList_total[[protein]]$H3N2_abFC,
                   plotList_total[[protein]]$Bvictoria_abFC,
                   plotList_total[[protein]]$Byamagata_abFC,
                   nrow = 1)
cowplot::plot_grid(plotList_total2[[protein]]$H1N1_abFC, 
                   plotList_total2[[protein]]$H3N2_abFC,
                   plotList_total2[[protein]]$Bvictoria_abFC,
                   plotList_total2[[protein]]$Byamagata_abFC,
                   nrow = 4)

# make pdf files
plotList_total3 <- list()
plotList_total4 <- list()
for (selectedVal in check_protein$OlinkID) {
  plotList3 <- list()
  plotList4 <- list()
  for(abType in c("H1N1_abFC", "H3N2_abFC", "Bvictoria_abFC", "Byamagata_abFC")) {
    plot_dat <- inputDat %>% filter(time == "T1") %>% 
      select(patientID, condition, disease, category, time, 
             all_of(c(abType, selectedVal))) %>%
      rename("abFoldchange" = 6, "expression" = 7)
    plotList3[[abType]] <- plot_dat %>% 
      ggscatter(x = "abFoldchange", y = "expression", 
                color = "disease", add = "reg.line") +
      stat_cor(aes(color = disease), method = "pearson") +
      xlab(paste0("log2_", abType)) + 
      ylab(paste0("log2_expression_", 
                  ZirFlu$protein_annot$Assay[which(ZirFlu$protein_annot$OlinkID == selectedVal) ]))
    
    plotList4[[abType]] <- plot_dat %>% 
      ggscatter(x = "abFoldchange", y = "expression", add = "reg.line") +
      stat_cor(method = "pearson")+
      xlab(paste0("log2_", abType)) + ylab(paste0("log2_expression_", 
                                                  ZirFlu$protein_annot$Assay[which(ZirFlu$protein_annot$OlinkID == selectedVal) ]))
  }
  plotList_total3[[selectedVal]] <- plotList3
  plotList_total4[[selectedVal]] <- plotList4
}

pdf(file = "protein_associAbFC_T1_v2.pdf", width = 12, onefile = TRUE)
for (protein in check_protein$OlinkID) {
  plot(cowplot::plot_grid(
    
    cowplot::plot_grid(plotList_total3[[protein]]$H1N1_abFC, 
                       plotList_total3[[protein]]$H3N2_abFC,
                       plotList_total3[[protein]]$Bvictoria_abFC,
                       plotList_total3[[protein]]$Byamagata_abFC,
                       nrow = 1),
    cowplot::plot_grid(plotList_total4[[protein]]$H1N1_abFC, 
                       plotList_total4[[protein]]$H3N2_abFC,
                       plotList_total4[[protein]]$Bvictoria_abFC,
                       plotList_total4[[protein]]$Byamagata_abFC,
                       nrow = 1),
    nrow = 2
  ))
}
dev.off()
# check if change to lm(abFC ~age + sex + proteinVal + disease) model 
protein <-  "OID20481" # IL17D 
protein <- "OID20575" # HSD11B1
protein <- "OID20624" # TNFSF12
cowplot::plot_grid(plotList_total3[[protein]]$H1N1_abFC, 
                   plotList_total3[[protein]]$H3N2_abFC,
                   plotList_total3[[protein]]$Bvictoria_abFC,
                   plotList_total3[[protein]]$Byamagata_abFC,
                   nrow = 2)
cowplot::plot_grid(plotList_total4[[protein]]$H1N1_abFC, 
                   plotList_total4[[protein]]$H3N2_abFC,
                   plotList_total4[[protein]]$Bvictoria_abFC,
                   plotList_total4[[protein]]$Byamagata_abFC,
                   nrow = 2)

# select consitent protein across season and strain
protein <- "OID20722" # FIS1
protein <- "OID20727" # GAL
protein <- "OID20596" # KRT19
protein <- "OID20669" # KYNU
protein <- "OID20768" # LTBR
protein <- "OID20480" # TNFRSF13C
protein <- "OID20611" # TNSFS10
protein <- "OID20533" # TRIM21
cowplot::plot_grid(plotList_total3[[protein]]$H1N1_abFC, 
                   plotList_total3[[protein]]$H3N2_abFC,
                   plotList_total3[[protein]]$Bvictoria_abFC,
                   plotList_total3[[protein]]$Byamagata_abFC,
                   nrow = 2)
cowplot::plot_grid(plotList_total4[[protein]]$H1N1_abFC, 
                   plotList_total4[[protein]]$H3N2_abFC,
                   plotList_total4[[protein]]$Bvictoria_abFC,
                   plotList_total4[[protein]]$Byamagata_abFC,
                   nrow = 2)

# check with the iMED data
load("/vol/projects/CIIM/Influenza/ZirrFlu/iMED_NhanNguyen/iMED.RData")

iMED_protein <- iMED$HAItiter %>%
  full_join(iMED$donorSample %>% 
  full_join(iMED$protein_dat %>% rownames_to_column("name")) %>%
  filter(time == "T1"))

check_protein_iMED <- iMED$protein_annot %>% 
  filter(Assay %in% c("IL17C", "IL17D", "TNFSF12", "TNFSF10", "IL6", 
                      "CXCL8", "SPRY2", "CD48", "HGF"))

iMED_plotList_total3 <- list()
iMED_plotList_total4 <- list()
for (selectedVal in check_protein_iMED$OlinkID) {
  plotList3 <- list()
  plotList4 <- list()
  for(abType in c("ab_H1N1", "ab_H3N2", "ab_B")) {
    plot_dat <- iMED_protein %>% 
      select(patientID, responder, all_of(c(abType, selectedVal))) %>%
      rename("abFoldchange" = 3, "expression" = 4)

    plotList3[[abType]] <- plot_dat %>% 
      ggscatter(x = "abFoldchange", y = "expression", 
                color = "responder", add = "reg.line") +
      stat_cor(aes(color = responder), method = "pearson") +
      xlab(paste0("log2_", abType)) + 
      ylab(paste0("log2_expression_", 
                  check_protein_iMED$Assay[which(check_protein_iMED$OlinkID == selectedVal) ]))
    
    plotList4[[abType]] <- plot_dat %>% 
      ggscatter(x = "abFoldchange", y = "expression", add = "reg.line") +
      stat_cor(method = "pearson") +
      xlab(paste0("log2_", abType)) + 
      ylab(paste0("log2_expression_", 
                  check_protein_iMED$Assay[which(check_protein_iMED$OlinkID == selectedVal) ]))
  }
  iMED_plotList_total3[[selectedVal]] <- plotList3
  iMED_plotList_total4[[selectedVal]] <- plotList4
}

pdf(file = "protein_associAbFC_T1_iMED.pdf", width = 12, onefile = TRUE)
for (protein in check_protein_iMED$OlinkID) {
  plot(
    cowplot::plot_grid(
      cowplot::plot_grid(iMED_plotList_total4[[protein]]$ab_H1N1, 
                         iMED_plotList_total4[[protein]]$ab_H3N2,
                         iMED_plotList_total4[[protein]]$ab_B, nrow = 1),
      cowplot::plot_grid(iMED_plotList_total3[[protein]]$ab_H1N1, 
                         iMED_plotList_total3[[protein]]$ab_H3N2,
                         iMED_plotList_total3[[protein]]$ab_B, nrow = 1),
      nrow = 2
    )
  )
}
dev.off()

protein <- "OID20722" # FIS1
protein <- "OID20727" # GAL
protein <- "OID20596" # KRT19
protein <- "OID20669" # KYNU
protein <- "OID20768" # LTBR
protein <- "OID20480" # TNFRSF13C
protein <- "OID20611" # TNSFS10
protein <- "OID20533" # TRIM21
cowplot::plot_grid(
  cowplot::plot_grid(iMED_plotList_total4[[protein]]$ab_H1N1, 
                     iMED_plotList_total4[[protein]]$ab_H3N2,
                     iMED_plotList_total4[[protein]]$ab_B, nrow = 1),
  cowplot::plot_grid(iMED_plotList_total3[[protein]]$ab_H1N1, 
                     iMED_plotList_total3[[protein]]$ab_H3N2,
                     iMED_plotList_total3[[protein]]$ab_B, nrow = 1),
  nrow = 2
)

# BOXPLOT OF EACH PROTEIN -------------
picked_proteins <- ZirFlu$protein_annot %>% 
  filter(Assay %in% c("TNFSF10", "CXCL8", "IL17D", "IL17C", "IL6"))
plot_dat <- inputDat %>% filter(time == "T1") %>% 
  select(patientID, condition, disease, category,
         picked_proteins$OlinkID)
datPlot <- plot_dat %>% 
  pivot_longer(cols = starts_with("OID"), names_to = "protein", values_to = "log2_proteinExp") %>%
  mutate(protein = gsub("OID20477", "IL17C", protein),
         protein = gsub("OID20481", "IL17D", protein),
         protein = gsub("OID20563", "IL6", protein),
         protein = gsub("OID20611", "TNFSF10", protein),
         protein = gsub("OID20624", "TNFSF12", protein),
         protein = gsub("OID20631", "CXCL8", protein)) # %>% filter(protein %in% c("IL6", "IL17C", "IL17D"))

bxp <- datPlot %>%
  ggboxplot(x = "protein", y = "log2_proteinExp", fill = "disease",
            bxp.errorbar = TRUE, palette = "npg",
            outlier.shape = 1, outlier.size = 4, outlier.color = "grey")
stat.test_condition <- datPlot %>% group_by(protein) %>%
  t_test(log2_proteinExp ~ disease)%>% add_xy_position(x = "protein", dodge = 0.8)

bxp + 
  stat_pvalue_manual(
    stat.test_condition, label = "p", tip.length = 0.01, bracket.nudge.y = 1)   + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) + ylim(0, 8.5)

load("/vol/projects/CIIM/Influenza/ZirrFlu/iMED_NhanNguyen/iMED.RData")
View(iMED$HAItiter)
View(iMED$protein_dat)
View(iMED$donorInfo)
View(iMED$donorSample)

dat <- iMED$donorSample %>% filter(time =="T1") %>%
  full_join(iMED$donorInfo) %>% 
  left_join(iMED$protein_dat %>% rownames_to_column("name"))

plotDat <- dat %>% select(patientID, responder, picked_proteins$OlinkID)

datPlot <- plotDat %>% 
  pivot_longer(cols = starts_with("OID"), names_to = "protein", values_to = "log2_proteinExp") %>%
  mutate(protein = gsub("OID20477", "IL17C", protein),
         protein = gsub("OID20481", "IL17D", protein),
         protein = gsub("OID20563", "IL6", protein),
         protein = gsub("OID20611", "TNFSF10", protein),
         protein = gsub("OID20624", "TNFSF12", protein),
         protein = gsub("OID20631", "CXCL8", protein)) %>%
  filter(protein %in% c("IL6", "IL17C", "IL17D")) %>%
  filter(responder %in% c("NR", "TR"))

bxp <- datPlot %>%
  ggboxplot(x = "protein", y = "log2_proteinExp", fill = "responder",
            bxp.errorbar = TRUE, palette = "npg",
            outlier.shape = 1, outlier.size = 4, outlier.color = "grey")
stat.test_condition <- datPlot %>% group_by(protein) %>%
  t_test(log2_proteinExp ~ responder)%>% add_xy_position(x = "protein", dodge = 0.8)

bxp + 
  stat_pvalue_manual(
    stat.test_condition, label = "p", tip.length = 0.01, bracket.nudge.y = 1)   + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) + ylim(0, 8.5)

# metabolite ----
inputDat <- ZirFlu$donorInfo %>% full_join(ZirFlu$HAItiter) %>%
  left_join(ZirFlu$donorSamples %>% filter(time == "T1")) %>%
  left_join(ZirFlu$metabolite_dat %>% rownames_to_column("probenID"))

check_metabolites <- ZirFlu$metabolite_annot %>% 
  filter(Formula %in% c("C20H30O5", "C20H34O5")) %>% 
  select(ionIdx, Formula) %>% distinct()

plotList_total3 <- list()
plotList_total4 <- list()
for (selectedVal in as.character(check_metabolites$ionIdx)) {
  plotList3 <- list()
  plotList4 <- list()
  for(abType in c("H1N1_abFC", "H3N2_abFC", "Bvictoria_abFC", "Byamagata_abFC")) {
    plot_dat <- inputDat %>% filter(time == "T1") %>% 
      select(patientID, condition, disease, category, time, 
             all_of(c(abType, selectedVal))) %>%
      rename("abFoldchange" = 6, "expression" = 7)
    plotList3[[abType]] <- plot_dat %>% 
      ggscatter(x = "abFoldchange", y = "expression", 
                color = "disease", add = "reg.line") +
      stat_cor(aes(color = disease), method = "pearson") +
      xlab(paste0("log2_", abType)) + 
      ylab(paste0("log2_expression_", 
                  check_metabolites$Formula[which(check_metabolites$ionIdx == selectedVal) ]))
    
    plotList4[[abType]] <- plot_dat %>% 
      ggscatter(x = "abFoldchange", y = "expression", add = "reg.line") +
      stat_cor(method = "pearson")+
      xlab(paste0("log2_", abType)) + 
      ylab(paste0("log2_expression_", 
                  check_metabolites$Formula[which(check_metabolites$ionIdx == selectedVal) ]))
  }
  plotList_total3[[selectedVal]] <- plotList3
  plotList_total4[[selectedVal]] <- plotList4
}

pdf(file = "metabolite_associAbFC_T1_v2.pdf", width = 12, onefile = TRUE)
for (protein in as.character(check_metabolites$ionIdx)) {
  plot(cowplot::plot_grid(
    
    cowplot::plot_grid(plotList_total3[[protein]]$H1N1_abFC, 
                       plotList_total3[[protein]]$H3N2_abFC,
                       plotList_total3[[protein]]$Bvictoria_abFC,
                       plotList_total3[[protein]]$Byamagata_abFC,
                       nrow = 1),
    cowplot::plot_grid(plotList_total4[[protein]]$H1N1_abFC, 
                       plotList_total4[[protein]]$H3N2_abFC,
                       plotList_total4[[protein]]$Bvictoria_abFC,
                       plotList_total4[[protein]]$Byamagata_abFC,
                       nrow = 1),
    nrow = 2
  ))
}
dev.off()

# check protein express overtime ------------------------------
inputDat <- ZirFlu$donorInfo %>% full_join(ZirFlu$donorSamples) %>%
  filter(patientID %in% ZirFlu$HAItiter$patientID) %>%
  left_join(ZirFlu$protein_dat %>% rownames_to_column("probenID")) %>%
  filter(season == "2019")

check_protein <- ZirFlu$protein_annot %>% 
  filter(OlinkID %in% selected_proteins)

plotList_total <- list()
for (protein in check_protein$OlinkID) {
  plotList <- list()
  
  plot_dat <- inputDat %>%
    select(patientID, condition, disease, time, {{protein}}) %>%
    pivot_wider(names_from = "time", values_from = protein) 
  
  plotList[["T2_vsT1"]] <- plot_dat %>% 
    ggscatter(x = "T1", y = "T2", 
              color = "disease", add = "reg.line") +
    stat_cor(aes(color = disease), method = "pearson") +
    xlab(paste0("log2Exp_T1_", 
                check_protein$Assay[which(check_protein$OlinkID == protein) ])) + 
    ylab(paste0("log2Exp_T2_",  
                check_protein$Assay[which(check_protein$OlinkID == protein) ]))
  
  plotList[["T3_vsT1"]] <- plot_dat %>% 
    ggscatter(x = "T1", y = "T3", 
              color = "disease", add = "reg.line") +
    stat_cor(aes(color = disease), method = "pearson") +
    xlab(paste0("log2Exp_T1_", 
                check_protein$Assay[which(check_protein$OlinkID == protein) ])) + 
    ylab(paste0("log2Exp_T3_", 
                check_protein$Assay[which(check_protein$OlinkID == protein) ]))
  
  plotList[["T3_vsT2"]] <- plot_dat %>% 
    ggscatter(x = "T2", y = "T3", 
              color = "disease", add = "reg.line") +
    stat_cor(aes(color = disease), method = "pearson") +
    xlab(paste0("log2Exp_T2_",  
                check_protein$Assay[which(check_protein$OlinkID == protein) ])) + 
    ylab(paste0("log2Exp_T3_",  
                check_protein$Assay[which(check_protein$OlinkID == protein) ]))
  
  plotList_total[[protein]] <- plotList
  
}

pdf(file = "overlapped_sigPro_overTime_2019.pdf", height = 4,width = 12, onefile = TRUE)
for (protein in check_protein$OlinkID) {
  plot(cowplot::plot_grid(plotList_total[[protein]]$T2_vsT1, 
                     plotList_total[[protein]]$T3_vsT1,
                     plotList_total[[protein]]$T3_vsT2,
                     nrow = 1))
}
dev.off()


# protein bet timepoint comparison ---

check_protein <- ZirFlu$protein_annot %>% 
  filter(OlinkID %in%  rownames(resPrePost_protein$`2019`$T1vsT3$padj))

check_protein <- ZirFlu$protein_annot %>% 
  filter(OlinkID %in%  rownames(resPrePost_protein$`2019`$T1vsT3$pval))

check_protein <- ZirFlu$protein_annot %>%
  filter(Assay %in% c("TNFSF10"))
# plotList_T2vsT1 <- list()
plotList_T3vsT1 <- list()
# plotList_T3vsT2 <- list()
for (protein in check_protein$OlinkID) {
  plot_dat <- inputDat %>%
    select(patientID, condition, disease, time, all_of(protein)) %>%
    pivot_wider(names_from = "time", values_from = protein) 
  
  # plotList_T2vsT1[[protein]] <- plot_dat %>% 
  #   ggscatter(x = "T1", y = "T2", 
  #             color = "disease", add = "reg.line") +
  #   stat_cor(aes(color = disease), method = "pearson") +
  #   xlab(paste0("log2Exp_T1_", 
  #               check_protein$Assay[which(check_protein$OlinkID == protein) ])) + 
  #   ylab(paste0("log2Exp_T2_",  
  #               check_protein$Assay[which(check_protein$OlinkID == protein) ]))
  
  plotList_T3vsT1[[protein]] <- plot_dat %>% 
    ggscatter(x = "T1", y = "T3", 
              color = "disease", add = "reg.line") +
    stat_cor(aes(color = disease), method = "pearson") +
    xlab(paste0("log2Exp_T1_", 
                check_protein$Assay[which(check_protein$OlinkID == protein) ])) + 
    ylab(paste0("log2Exp_T3_", 
                check_protein$Assay[which(check_protein$OlinkID == protein) ]))
  
  # plotList_T3vsT2[[protein]] <- plot_dat %>% 
  #   ggscatter(x = "T2", y = "T3", 
  #             color = "disease", add = "reg.line") +
  #   stat_cor(aes(color = disease), method = "pearson") +
  #   xlab(paste0("log2Exp_T2_",  
  #               check_protein$Assay[which(check_protein$OlinkID == protein) ])) + 
  #   ylab(paste0("log2Exp_T3_",  
  #               check_protein$Assay[which(check_protein$OlinkID == protein) ]))
}

library(gridExtra)
# pdf(file = "selectPro_T2vsT1.pdf", width = 12, onefile = TRUE)
# ml1 <- marrangeGrob(plotList_T2vsT1, nrow = 2, ncol = 3)
# ml1
# dev.off()

# pdf(file = "selectPro_T3vsT1.pdf", width = 12, onefile = TRUE)
# ml1 <- marrangeGrob(plotList_T3vsT1, nrow = 2, ncol = 3)
# ml1
# dev.off()
# 
# pdf(file = "selectPro_T2vsT1.pdf", width = 12, onefile = TRUE)
# ml1 <- marrangeGrob(plotList_T2vsT1, nrow = 2, ncol = 3)
# ml1
# dev.off()

pdf(file = "selectPro_sigPadj_T3vsT1.pdf", width = 12, onefile = TRUE)
ml1 <- marrangeGrob(plotList_T3vsT1, nrow = 2, ncol = 3)
ml1
dev.off()

pdf(file = "selectPro_sigPval_T3vsT1.pdf", width = 12, onefile = TRUE)
ml1 <- marrangeGrob(plotList_T3vsT1, nrow = 2, ncol = 3)
ml1
dev.off()

ggboxplot(plot_dat, x = "time", y = "OID20611", 
          palette = "npg", add = "jitter", color = "disease")

check_protein <- ZirFlu$protein_annot %>% 
  filter(OlinkID %in%  rownames(resPrePost_protein$`2019`$T1vsT3$pval))
plotList <- list()
for (protein in check_protein$OlinkID) {
  plot_dat <- inputDat %>%
    select(patientID, condition, disease, time, all_of(protein)) %>%
    rename("log2_proExp" = protein)
  
  plotList[[protein]] <- ggboxplot(plot_dat, x = "time", y = "log2_proExp", 
                                   palette = "npg", add = "jitter", color = "disease") +
  ylab(paste0("log2_proExp_", check_protein$Assay[which(check_protein$OlinkID == protein) ]))
}

pdf(file = "selectPro_sigPval_T3vsT1_boxplot.pdf", width = 12, onefile = TRUE)
ml1 <- marrangeGrob(plotList, nrow = 2, ncol = 3)
ml1
dev.off()


# correlation betweel change in abTiter vs. change in protein
inputDat <- ZirFlu$donorInfo %>% full_join(ZirFlu$HAItiter) %>%
  left_join(ZirFlu$donorSamples %>% filter(time != "T2")) %>%
  left_join(ZirFlu$protein_dat %>% rownames_to_column("probenID")) %>%
  filter(season == "2019")

ab_T3vsT1 <- inputDat %>%
  select(patientID, condition, disease, matches("T1|T3")) %>%
  mutate(H1N1_T3vsT1 = H1N1_T3 - H1N1_T1) %>%
  mutate(H3N2_T3vsT1 = H3N2_T3 - H3N2_T1) %>%
  mutate(Bvictoria_T3vsT1 = Bvictoria_T3 - Bvictoria_T1) %>%
  mutate(Byamagata_T3vsT1 = Byamagata_T3 - Byamagata_T1) %>% distinct()

plot_dat <- inputDat %>%
  select(patientID, condition, disease, time, all_of(protein)) %>%
  pivot_wider(names_from = "time", values_from = protein) %>%
  mutate(protein_T3vsT1 = T3 - T1)

plot_dat2 <- ab_T3vsT1 %>% full_join(plot_dat)

plot_dat2 %>% 
  ggscatter(x = "H1N1_T3vsT1", y = "protein_T3vsT1", 
            color = "disease", add = "reg.line") +
  stat_cor(aes(color = disease), method = "pearson") +
  xlab("log2_H1N1_T3vsT1") + 
  ylab(paste0("log2Exp_proT3vsT1_", 
              check_protein$Assay[which(check_protein$OlinkID == protein) ]))

plotList_total3 <- list()
plotList_total4 <- list()
for (protein in check_protein$OlinkID) {
  plotList3 <- list()
  plotList4 <- list()
  
  plot_dat <- inputDat %>%
    select(patientID, condition, disease, time, all_of(protein)) %>%
    pivot_wider(names_from = "time", values_from = protein) %>%
    mutate(protein_T3vsT1 = T3 - T1)
  
  plot_dat2 <- ab_T3vsT1 %>% full_join(plot_dat)
  
  for(abTiter in c("H1N1_T3vsT1", "H3N2_T3vsT1", "Bvictoria_T3vsT1", "Byamagata_T3vsT1")) {

    
    plotList3[[abTiter]] <- plot_dat2 %>% 
      ggscatter(x = abTiter, y = "protein_T3vsT1", 
                color = "disease", add = "reg.line") +
      stat_cor(aes(color = disease), method = "pearson") +
      xlab("log2_H1N1_T3vsT1") + 
      ylab(paste0("log2Exp_proT3vsT1_", 
                  check_protein$Assay[which(check_protein$OlinkID == protein) ]))
    
    plotList4[[abTiter]] <- plot_dat2 %>% 
      ggscatter(x = abTiter, y = "protein_T3vsT1", add = "reg.line") +
      stat_cor(method = "pearson") +
      xlab("log2_H1N1_T3vsT1") + 
      ylab(paste0("log2Exp_proT3vsT1_", 
                  check_protein$Assay[which(check_protein$OlinkID == protein) ]))
  }
  plotList_total3[[protein]] <- plotList3
  plotList_total4[[protein]] <- plotList4
}

pdf(file = "proPval_T3vsT1_2019.pdf", width = 12, onefile = TRUE)
for (protein in check_protein$OlinkID) {
  plot(cowplot::plot_grid(
    
    cowplot::plot_grid(plotList_total3[[protein]]$H1N1_T3vsT1, 
                       plotList_total3[[protein]]$H3N2_T3vsT1,
                       plotList_total3[[protein]]$Bvictoria_T3vsT1,
                       plotList_total3[[protein]]$Byamagata_T3vsT1,
                       nrow = 1),
      cowplot::plot_grid(plotList_total4[[protein]]$H1N1_T3vsT1, 
                       plotList_total4[[protein]]$H3N2_T3vsT1,
                       plotList_total4[[protein]]$Bvictoria_T3vsT1,
                       plotList_total4[[protein]]$Byamagata_T3vsT1,
                       nrow = 1),
    nrow = 2
  ))
}
dev.off()

# check protein expression with liver disease severities -----------------------
inputDat <- inputDat %>% left_join(a) # check the code "clean_ValerieInfo"

# make pdf files
variables <- c("CHILD-Pugh ponts",  "CHILD-Pugh state of cirrhosis (A=1, B=2, C=3)",
               "MELD points", "esophageal_varices")
check_protein <- ZirFlu$protein_annot %>% 
  filter(Assay %in% c("TNFSF10", "CXCL8", "TNFSF12", "IL6", "IL17C", "IL17D", "FOXO1", "PLXNA4"))
plotList_total3 <- list()
plotList_total4 <- list()
for (selectedVal in check_protein$OlinkID) {
  plotList3 <- list()
  plotList4 <- list()
  for(variable in variables) {
    plot_dat <- inputDat %>% filter(time == "T1") %>% 
      select(patientID, condition, disease, category, time, 
             all_of(c(variable, selectedVal))) %>%
      rename("variable" = 6, "expression" = 7)
    plotList3[[variable]] <- plot_dat %>% 
      ggscatter(x = "variable", y = "expression", 
                color = "condition", add = "reg.line") +
      stat_cor(aes(color = condition), method = "pearson") +
      xlab(variable) + 
      ylab(paste0("log2_expression_", 
                  ZirFlu$protein_annot$Assay[which(ZirFlu$protein_annot$OlinkID == selectedVal) ]))
    
    plotList4[[variable]] <- plot_dat %>% 
      ggscatter(x = "variable", y = "expression", add = "reg.line") +
      stat_cor(method = "pearson")+
      xlab(variable) + ylab(paste0("log2_expression_", 
                                                  ZirFlu$protein_annot$Assay[which(ZirFlu$protein_annot$OlinkID == selectedVal) ]))
  }
  plotList_total3[[selectedVal]] <- plotList3
  plotList_total4[[selectedVal]] <- plotList4
}

pdf(file = "protein_associCirrhosisScore_T1.pdf", width = 12, onefile = TRUE)
for (protein in check_protein$OlinkID) {
  plot(cowplot::plot_grid(
    
    cowplot::plot_grid(plotList_total3[[protein]]$"CHILD-Pugh ponts", 
                       plotList_total3[[protein]]$"CHILD-Pugh state of cirrhosis (A=1, B=2, C=3)",
                       plotList_total3[[protein]]$"MELD points",
                       plotList_total3[[protein]]$"esophageal_varices",
                       nrow = 1),
    cowplot::plot_grid(plotList_total4[[protein]]$"CHILD-Pugh ponts", 
                       plotList_total4[[protein]]$"CHILD-Pugh state of cirrhosis (A=1, B=2, C=3)",
                       plotList_total4[[protein]]$"MELD points",
                       plotList_total4[[protein]]$"esophageal_varices",
                       nrow = 1),
    nrow = 2
  ))
}
dev.off()
