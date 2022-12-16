library(tidyverse)
library(lme4)
library(lmerTest)
library(broom.mixed)
library(ggpubr)
library(ggVennDiagram)
library(ComplexHeatmap)

# note switich the lm() model from proteinVal ~sex + age + abTiterFC + disease 
# to abTiterFC ~sex + age + proteinValue + disease
# specific function for this code ---------------------------
get.lm_forDat_v2 <- function(inputDat, varList) {
  lm_results <- list()
  abTiters <- c("H1N1_abFC", "H3N2_abFC", 
                "Bvictoria_abFC", "Byamagata_abFC")
  for (name in abTiters) {
    lm_results[[name]] <- get.lmTest_correction_v2(abFC = name,
                                                   variableList = varList, inputDat = inputDat)
  }
  return(lm_results)
}

get.lm_padj <- function(lm_results) {
  lm_padj <- list()
  for (i in names(lm_results)) {
    dat_temp <- p.adjust(lm_results[[i]]$p.value, method = "fdr")
    lm_padj[[i]] <- cbind(lm_results[[i]], "p.adj" = dat_temp)
  }
  return(lm_padj)
}

get.protein_associ_abFC <- function(lm_results) {
  protein_associ_abFC <- list()
  
  for (i in names(lm_results)) {
    protein_associ_abFC[[i]] <- lm_results[[i]] %>% filter(p.value < 0.05) %>%
      filter(independentVariable == "proteinVal")
  }
  return(protein_associ_abFC)
}

get.protein_associ_abFC_adj <- function(lm_results) {
  protein_associ_abFC_adj <- list()
  
  for (i in names(lm_results)) {
    protein_associ_abFC_adj[[i]] <- lm_results[[i]] %>% filter(p.adj < 0.05) %>%
      filter(independentVariable == "proteinVal")
  }
  return(protein_associ_abFC_adj)
}

# protein with lm() model ---------------------------------------------------
inputDat_protein <- ZirFlu$donorSamples %>% 
  full_join(ZirFlu$donorInfo) %>% filter(time == "T1") %>%
  filter(patientID %in% ZirFlu$HAItiter$patientID) %>%
  left_join(ZirFlu$protein_dat %>% rownames_to_column("probenID"))

# per season 
input_protein <- list()
input_protein[["2019"]] <- ZirFlu$HAItiter_2019 %>% 
  left_join(inputDat_protein %>% filter(season == "2019"))
input_protein[["2020"]] <- ZirFlu$HAItiter_2020 %>% 
  left_join(inputDat_protein %>% filter(season == "2020"))

lmRes_protein <- list()
for (season in c("2019", "2020")) {
  lmRes_protein[[season]]$lmRes <- get.lm_forDat_v2(input_protein[[season]], 
                                                 varList = names(ZirFlu$protein_dat))
  lmRes_protein[[season]]$lmAdj <- lmRes_protein[[season]]$lmRes %>% get.lm_padj()
  lmRes_protein[[season]]$sig <- lmRes_protein[[season]]$lmRes %>% get.protein_associ_abFC()
  lmRes_protein[[season]]$sigAdj <- lmRes_protein[[season]]$lmAdj %>% get.protein_associ_abFC_adj()
}

# save result
save(lmRes_protein, file = "temp/lmResult_protein_2seasons_v2.RData")

# check protein
get.sigPro_acrossSeason <- function(lmRes, strain) {
  outcome <- intersect(lmRes$`2019`$sig[[strain]]$targetVariable,
                       lmRes$`2020`$sig[[strain]]$targetVariable)
  return(outcome)
}

proSig <- list()
for (strain in c("H1N1_abFC", "H3N2_abFC", "Bvictoria_abFC", "Byamagata_abFC")) {
  proSig[[strain]] <- get.sigPro_acrossSeason(lmRes_protein, strain)
}

unlist(proSig)
get.proteinAnnot(ZirFlu$protein_annot, unlist(proSig))

selected_proteins <- unique(c(lmRes_protein$`2019`$sig$H1N1_abFC$targetVariable,
       lmRes_protein$`2019`$sig$H3N2_abFC$targetVariable,
       lmRes_protein$`2019`$sig$Bvictoria_abFC$targetVariable,
       lmRes_protein$`2019`$sig$Byamagata_abFC$targetVariable,
       lmRes_protein$`2020`$sig$H1N1_abFC$targetVariable,
       lmRes_protein$`2020`$sig$H3N2_abFC$targetVariable,
       lmRes_protein$`2020`$sig$Bvictoria_abFC$targetVariable,
       lmRes_protein$`2020`$sig$Byamagata_abFC$targetVariable))


# make heatmap
get.lmStatistic <- function(inputDat) {
  outcome <- inputDat %>% 
    lapply(function(x) x %>% select(targetVariable, statistic)) %>% 
    bind_rows(.id = "groups") %>% mutate(groups = substring(groups, 1, 7)) %>%
    pivot_wider(names_from = groups, values_from = statistic) %>% 
    column_to_rownames("targetVariable")
  return(outcome)
}

plotDat_wide <- list()
plotDat_long <- list()
for (season in c("2019", "2020")) {
  varlmAbTiter <- lmRes_protein[[season]]$lmRes %>% 
    lapply(function(x) x%>% filter(independentVariable == "proteinVal")) %>%
    get.lmStatistic()
  
  heatmapDat <- varlmAbTiter[which(rownames(varlmAbTiter) %in% selected_proteins),]
  
  plotDat_wide[[season]] <- heatmapDat %>% #as.data.frame %>% 
    rownames_to_column("OlinkID") %>% left_join(ZirFlu$protein_annot) %>%
    select(-c(OlinkID, UniProt)) %>% column_to_rownames("Assay")
  
  plotDat_long[[season]] <- heatmapDat %>% #as.data.frame() %>% 
    rownames_to_column("OlinkID") %>% left_join(ZirFlu$protein_annot) %>%
    pivot_longer(cols = c(2:5), names_to = "time", values_to = "statistic") #%>%
    # get.dat_sigPvalue(., lm_sigPvalue = lmRes_protein[[season]]$sig$var_onlyAbTiter) %>%
    # rename("significance_pvalue" = "significance") %>% 
    # get.dat_sigPvalue(., lm_sigPvalue = lmRes_protein[[season]]$sigAdj$var_onlyAbTiter)
}

plotDat_wide_2seasons <- plotDat_wide$`2019` %>% rownames_to_column("protein") %>%
  rename("2019_H1N1" = "H1N1_ab", "2019_H3N2" = "H3N2_ab",
         "2019_Bvic" = "Bvictor", "2019_Byam"= "Byamaga") %>% 
  full_join(plotDat_wide$`2020` %>% rownames_to_column("protein") %>%
              rename("2020_H1N1" = "H1N1_ab", "2020_H3N2" = "H3N2_ab",
                     "2020_Bvic" = "Bvictor", "2020_Byam" = "Byamaga")) %>%
  column_to_rownames("protein")

plotDat_wide_2seasons %>% as.matrix() %>% heatmap(Colv = NA, Rowv = NA)
plotDat_wide_2seasons %>% as.matrix() %>% Heatmap()
plotDat_wide_2seasons %>% as.matrix() %>% Heatmap(cluster_columns = FALSE)

plotDat_long_2seasons <- plotDat_long$`2019` %>% mutate(season = "2019") %>%
  rbind(plotDat_long$`2020` %>% mutate(season = "2020")) %>%
  unite(time, season, time) %>% mutate(time = substring(time, 1, 9))

plotDat_long_2seasons  %>% 
  ggplot(aes(x = time, y = Assay)) + 
  geom_tile(aes(fill = statistic)) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plotDat_long_2seasons  %>% 
  ggplot(aes(x = time, y = Assay)) + 
  geom_tile(aes(fill = statistic)) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plotDat_long_2seasons %>% 
  ggplot(aes(x = time, y = Assay)) + 
  geom_tile(aes(fill = statistic)) +
  scale_y_discrete(limits = plotDat_long_2seasons$Assay[hclust(dist(plotDat_wide_2seasons))$order])+
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme_bw() +
  scale_y_discrete(labels = labels)

# new heat map with consistent within strain  across season
plotDat_long_2seasons$time <- factor(plotDat_long_2seasons$time,
                                     c("2019_H1N1", "2020_H1N1",
                                       "2019_H3N2", "2020_H3N2",
                                       "2019_Bvic", "2020_Bvic",
                                       "2019_Byam", "2020_Byam"))
highlight_proteins <- plotDat_long_2seasons %>% group_by(Assay) %>%
  summarise(postitive = sum(statistic > 0, na.rm = TRUE)) %>%
  full_join(plotDat_long_2seasons %>% group_by(Assay) %>%
              summarise(negative = sum(statistic < 0, na.rm = TRUE))) %>%
  mutate(select = ifelse(postitive >= 8*0.75 | negative >= 8*0.75, TRUE, NA)) %>%
  filter(select == TRUE)

highlight_proteins <- plotDat_long_2seasons %>% group_by(Assay) %>%
  summarise(postitive = sum(statistic > 0, na.rm = TRUE)) %>%
  full_join(plotDat_long_2seasons %>% group_by(Assay) %>%
              summarise(negative = sum(statistic < 0, na.rm = TRUE))) %>%
  mutate(select = ifelse(postitive >= 8 | negative >= 8, TRUE, NA)) %>%
  filter(select == TRUE)

plotDat_long_2seasons  %>% filter(Assay %in% highlight_proteins$Assay) %>%
  ggplot(aes(x = time, y = Assay)) + 
  geom_tile(aes(fill = statistic)) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

