library(tidyverse)
library(ggpubr)
library(rstatix)
# function
get.plot_condition_category <- function(dat, var, abTiter) {
  plot_dat <- dat %>% 
    rename("var" = all_of(var), "abTiter" = all_of(abTiter))
  
  ggplot(plot_dat, aes(x = abTiter, y = var)) +
    geom_point(aes(color = category, shape = condition), 
               size = 3, alpha = 0.8,
               position = position_jitter(w = 0.08, h = 0.1)) +
    geom_smooth(method = "lm", se = FALSE) + theme_bw() +
    xlab(paste0("log2_", abTiter)) + ylab(paste0("log2_", var))
}

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
# input data  ------------------------------------------------------------------
inputDat <- ZirFlu$donorInfo %>% full_join(ZirFlu$HAItiter) %>%
  left_join(ZirFlu$donorSamples %>% filter(time == "T1")) %>%
  left_join(ZirFlu$protein_dat %>% rownames_to_column("probenID")) #%>% filter(season == "2020")

inputDat$category <- factor(inputDat$category, levels = c("NR", "Other", "TR"))
# T1 protein plot -----------------------------------------------------------
## scater plot with lm()
ggplot(inputDat, aes(x = H1N1_abFC, y = OID20481)) +
  geom_point(aes(color = category, shape = condition), 
             size = 3, alpha = 0.8,
             position = position_jitter(w = 0.08, h = 0.1)) +
  geom_smooth(method = "lm", se = FALSE) + theme_bw()

get.plot_condition_category(inputDat, var = "OID20481", abTiter = "H1N1_abFC")

cowplot::plot_grid(
  get.plot_condition_category(inputDat, var = "OID20481", abTiter = "H1N1_abFC"),
  get.plot_condition_category(inputDat, var = "OID20624", abTiter = "H1N1_abFC"),
  get.plot_condition_category(inputDat, var = "OID20765", abTiter = "H1N1_abFC"),
  nrow = 1
)

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
  nrow = 1
)
## boxplot
cowplot::plot_grid(
  ggboxplot(inputDat, x = "category", y = "OID20481", palette = "jco", add = "jitter"),
  ggboxplot(inputDat, x = "disease", y = "OID20481", palette = "jco", add = "jitter"),
  ggboxplot(inputDat, x = "condition", y = "OID20481", palette = "jco", add = "jitter") + 
    rotate_x_text(20),
  ggboxplot(inputDat, x = "category", y = "OID20624", palette = "jco", add = "jitter"),
  ggboxplot(inputDat, x = "disease", y = "OID20624", palette = "jco", add = "jitter"),
  ggboxplot(inputDat, x = "condition", y = "OID20624", palette = "jco", add = "jitter") + 
    rotate_x_text(20),
  ggboxplot(inputDat, x = "category", y = "OID20765", palette = "jco", add = "jitter"),
  ggboxplot(inputDat, x = "disease", y = "OID20765", palette = "jco", add = "jitter"),
  ggboxplot(inputDat, x = "condition", y = "OID20765", palette = "jco", add = "jitter") + 
    rotate_x_text(20),
  nrow = 3, byrow = FALSE
)

## check multiple proteins
#check_protein <- ZirFlu$protein_annot %>% filter(Olink %in% selected_proteins)
selected_proteins <- c("TNFSF10", "CXCL8", "IL17D", "IL17D", "IL6")
check_protein <- ZirFlu$protein_annot %>% filter(Assay %in% selected_proteins)

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

# check proteins
protein <- "OID20611" # TNFSF10
cowplot::plot_grid(plotList_total[[protein]]$H1N1_abFC, 
                   plotList_total[[protein]]$H3N2_abFC,
                   plotList_total[[protein]]$Bvictoria_abFC,
                   plotList_total[[protein]]$Byamagata_abFC,
                   nrow = 1)
cowplot::plot_grid(plotList_total2[[protein]]$H1N1_abFC, 
                   plotList_total2[[protein]]$H3N2_abFC,
                   plotList_total2[[protein]]$Bvictoria_abFC,
                   plotList_total2[[protein]]$Byamagata_abFC,
                   nrow = 1)
cowplot::plot_grid(plotList_total3[[protein]]$H1N1_abFC, 
                   plotList_total3[[protein]]$H3N2_abFC,
                   plotList_total3[[protein]]$Bvictoria_abFC,
                   plotList_total3[[protein]]$Byamagata_abFC,
                   nrow = 1)
cowplot::plot_grid(plotList_total4[[protein]]$H1N1_abFC, 
                   plotList_total4[[protein]]$H3N2_abFC,
                   plotList_total4[[protein]]$Bvictoria_abFC,
                   plotList_total4[[protein]]$Byamagata_abFC,
                   nrow = 1)

## make the plot for the manuscript TFNSF10 & CXCL8
protein <- "OID20611" # TNFSF10
protein <- "OID20631" # CXCL8
cowplot::plot_grid(
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
)

## make pdf files of the protein expression -------------------------------------
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
## check with the iMED data ----------------------------------------------------
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

# metabolite --------------------------------------------
inputDat <- ZirFlu$donorInfo %>% full_join(ZirFlu$HAItiter) %>%
  left_join(ZirFlu$donorSamples %>% filter(time == "T1")) %>%
  left_join(ZirFlu$metabolite_dat %>% rownames_to_column("probenID"))# %>% filter(season == "2019")

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

# box plot for each metabolites  -------------
picked_metabolites <- ZirFlu$metabolite_annot %>% select(ionIdx, Formula) %>%
  filter(ionIdx %in% c("782", "935", "944")) %>% distinct()
plot_dat <- inputDat %>% 
  select(patientID, season, time, condition, disease, as.character(picked_metabolites$ionIdx)) %>%
  rename("C20H32O2" = "782", "C20H30O5" = "935", "C20H34O5" = "944")

datPlot <- plot_dat %>% 
  pivot_longer(cols = starts_with("C20"), names_to = "metabolite", values_to = "log2_metaExp")

bxp <- datPlot %>%
  ggboxplot(x = "metabolite", y = "log2_metaExp", fill = "disease",
            bxp.errorbar = TRUE, palette = "npg",
            outlier.shape = 1, outlier.size = 4, outlier.color = "grey")
stat.test_condition <- datPlot %>% group_by(metabolite) %>%
  t_test(log2_metaExp ~ disease) %>% add_xy_position(x = "metabolite", dodge = 0.8)

bxp + 
  stat_pvalue_manual(
    stat.test_condition, label = "p", tip.length = 0.01, bracket.nudge.y = 1)   + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) + ylim(13, 20)

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

