library(tidyverse)
library(lme4)
library(lmerTest)
library(broom.mixed)
library(ggpubr)
library(ggVennDiagram)
library(ComplexHeatmap)
library("KEGGREST")
# specific function for this code ---------------------------
get.lm_forDat <- function(inputDat, varList) {
  lm_results <- list()
  abTiters <- c("H1N1_abFC", "H3N2_abFC", 
                "Bvictoria_abFC", "Byamagata_abFC")
  for (name in abTiters) {
    lm_results[[name]] <- get.lmTest_correction(variableList = varList,
                                                abTiter = name, inputDat = inputDat)
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

get.var_associFactor <- function(lm_results) {
  var_associDisease <- list()
  var_associAbTiter <- list()
  var_onlyAbTiter <- list()
  
  for (i in names(lm_results)) {
    var_associDisease[[i]] <- lm_results[[i]] %>% filter(p.value < 0.05) %>%
      filter(independentVariable == "diseasehealthy")
    
    var_associAbTiter[[i]] <- lm_results[[i]] %>% filter(p.value < 0.05) %>%
      filter(independentVariable == "abTiter")
    
    var_onlyAbTiter[[i]] <- get.var_onlyAbTiter(var_associAbTiter, 
                                                var_associDisease,
                                                abTiter = i)
  }
  return(list("var_associDisease" = var_associDisease, 
              "var_associAbTiter" = var_associAbTiter, 
              "var_onlyAbTiter" = var_onlyAbTiter))
}
get.var_associFactor_padj <- function(lm_results, lmRes_pvalue) {
  var_associDisease <- list()
  var_associAbTiter <- list()
  var_onlyAbTiter <- list()
  
  for (i in names(lm_results)) {
    var_associDisease[[i]] <- lm_results[[i]] %>% filter(p.adj < 0.05) %>%
      filter(independentVariable == "diseasehealthy")
    
    var_associAbTiter[[i]] <- lm_results[[i]] %>% filter(p.adj < 0.05) %>%
      filter(independentVariable == "abTiter")
    
    var_onlyAbTiter_temp <- get.var_onlyAbTiter(var_associAbTiter, 
                                                var_associDisease,
                                                abTiter = i)
    var_onlyAbTiter[[i]] <- var_onlyAbTiter_temp[which(var_onlyAbTiter_temp %in% lmRes_pvalue[[i]])]
  }
  return(list("var_associDisease" = var_associDisease, 
              "var_associAbTiter" = var_associAbTiter, 
              "var_onlyAbTiter" = var_onlyAbTiter))
}

get.targetedVar <- function(inputDat) {
  outcome <- inputDat %>% 
    lapply(function(x) x %>% select(targetVariable))
  return(outcome)
}

# make heatmap
get.lmStatistic <- function(inputDat) {
  outcome <- inputDat %>% 
    lapply(function(x) x %>% select(targetVariable, statistic)) %>% 
    bind_rows(.id = "groups") %>% mutate(groups = substring(groups, 1, 7)) %>%
    pivot_wider(names_from = groups, values_from = statistic) %>% 
    column_to_rownames("targetVariable")
  return(outcome)
}

get.dat_sigPvalue <- function(inputDat, lm_sigPvalue) {
  outcome <- as.data.frame(matrix(nrow=0, ncol = 4))
  for (name in names(lm_sigPvalue)) {
    dat_temp <- inputDat %>% filter(time == substring(name, 1, 7)) %>% 
      mutate(significance = ifelse(OlinkID %in% lm_sigPvalue[[name]]$targetVariable, TRUE, NA))
    outcome <- rbind(outcome, dat_temp)
  }
  return(outcome)
}

get.dat_sigPvalue_metabolite <- function(inputDat, lm_sigPvalue) {
  outcome <- as.data.frame(matrix(nrow=0, ncol = 4))
  for (name in names(lm_sigPvalue)) {
    dat_temp <- inputDat %>% filter(time == substring(name, 1, 7)) %>% 
      mutate(significance = ifelse(ionIdx %in% lm_sigPvalue[[name]], TRUE, NA))
    outcome <- rbind(outcome, dat_temp)
  }
  return(outcome)
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
  lmRes_protein[[season]]$lmRes <- get.lm_forDat(input_protein[[season]], 
                                                 varList = names(ZirFlu$protein_dat))
  lmRes_protein[[season]]$lmAdj <- lmRes_protein[[season]]$lmRes %>% get.lm_padj()
  lmRes_protein[[season]]$sig <- lmRes_protein[[season]]$lmRes %>% get.var_associFactor()
  lmRes_protein[[season]]$sigAdj <- lmRes_protein[[season]]$lmAdj %>% 
    get.var_associFactor_padj(lmRes_protein[[season]]$sig$var_onlyAbTiter)
}
# save result
save(lmRes_protein, file = "temp/lmResult_protein_2seasons.RData")

## downstream analysis  ---------------------------------------------------------
load("temp/lmResult_protein_2seasons.RData")

## check top protein ---
proSig_onlyAb_padj_2019 <- unique(unlist(lmRes_protein$`2019`$sigAdj$var_onlyAbTiter))
proSig_onlyAb_padj_2020 <- unique(unlist(lmRes_protein$`2020`$sigAdj$var_onlyAbTiter))

proSig_onlyAb_2019 <- unique(unlist(lmRes_protein$`2019`$sig$var_onlyAbTiter))
proSig_onlyAb_2020 <- unique(unlist(lmRes_protein$`2020`$sig$var_onlyAbTiter))

proSig_ab_2019 <- unique(c(lmRes_protein$`2019`$sig$var_associAbTiter$H1N1_abFC$targetVariable,
                           lmRes_protein$`2019`$sig$var_associAbTiter$H3N2_abFC$targetVariable,
                           lmRes_protein$`2019`$sig$var_associAbTiter$Bvictoria_abFC$targetVariable,
                           lmRes_protein$`2019`$sig$var_associAbTiter$Byamagata_abFC$targetVariable))
proSig_ab_2020 <- unique(c(lmRes_protein$`2020`$sig$var_associAbTiter$H1N1_abFC$targetVariable,
                           lmRes_protein$`2020`$sig$var_associAbTiter$H3N2_abFC$targetVariable,
                           lmRes_protein$`2020`$sig$var_associAbTiter$Bvictoria_abFC$targetVariable,
                           lmRes_protein$`2020`$sig$var_associAbTiter$Byamagata_abFC$targetVariable))

get.proteinAnnot(ZirFlu$protein_annot, c(proSig_onlyAb_padj_2019, proSig_onlyAb_padj_2020)) 
# heatmap --------------------
selected_proteins <- unique(c(proSig_onlyAb_padj_2019, proSig_onlyAb_padj_2020))
selected_proteins <- unique(c(proSig_onlyAb_2019, proSig_onlyAb_2020))

sigProtein_ab <- unique(c(proSig_ab_2019, proSig_ab_2020))
sigProtein_onlyAb <- unique(c(proSig_onlyAb_2019, proSig_onlyAb_2020))
selected_proteins <- sigProtein_ab[-which(sigProtein_ab %in% sigProtein_onlyAb)]

# picked_proteins <- ZirFlu$protein_annot %>% 
#   filter(Assay %in% c("TNFSF10", "TNFSF12", "CXCL8", "IL6", "IL17C", "IL17D"))
# selected_proteins <- picked_proteins$OlinkID

plotDat_wide <- list()
plotDat_long <- list()
for (season in c("2019", "2020")) {
  varlmAbTiter <- lmRes_protein[[season]]$lmRes %>% 
    lapply(function(x) x%>% filter(independentVariable == "abTiter")) %>%
    get.lmStatistic()
  
  heatmapDat <- varlmAbTiter[which(rownames(varlmAbTiter) %in% selected_proteins),]
  
  plotDat_wide[[season]] <- heatmapDat %>%
    rownames_to_column("OlinkID") %>% left_join(ZirFlu$protein_annot) %>%
    select(-c(OlinkID, UniProt)) %>% column_to_rownames("Assay")
  
  plotDat_long[[season]] <- heatmapDat %>%
    rownames_to_column("OlinkID") %>% left_join(ZirFlu$protein_annot) %>%
    pivot_longer(cols = c(2:5), names_to = "time", values_to = "statistic") %>%
    get.dat_sigPvalue(., lm_sigPvalue = lmRes_protein[[season]]$sig$var_associAbTiter) %>%
    rename("significance_pvalue" = "significance") %>% 
    get.dat_sigPvalue(., lm_sigPvalue = lmRes_protein[[season]]$sigAdj$var_associAbTiter)
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

# # if need heat map with consistent within strain  across season
# plotDat_long_2seasons$time <- factor(plotDat_long_2seasons$time,
#                                      c("2019_H1N1", "2020_H1N1",
#                                        "2019_H3N2", "2020_H3N2",
#                                        "2019_Bvic", "2020_Bvic",
#                                        "2019_Byam", "2020_Byam"))
plotDat_long_2seasons  %>% 
  ggplot(aes(x = time, y = Assay)) + 
  geom_tile(aes(fill = statistic)) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") + 
  geom_text(aes(label = ifelse(significance == TRUE, "**", ""))) +
  geom_text(aes(label = ifelse(significance_pvalue == TRUE, "*", ""))) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

plotDat_long_2seasons %>% 
  ggplot(aes(x = time, y = Assay)) + 
  geom_tile(aes(fill = statistic)) +
  scale_y_discrete(limits = plotDat_long_2seasons$Assay[hclust(dist(plotDat_wide_2seasons))$order])+
  geom_text(aes(label = ifelse(significance == TRUE, "**", ""))) +
  geom_text(aes(label = ifelse(significance_pvalue == TRUE, "*", ""))) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") + theme_bw() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

# only select protein show consistent trend (>75% across strains and years)
highlight_proteins <- plotDat_long_2seasons %>% group_by(Assay) %>%
  summarise(postitive = sum(statistic > 0, na.rm = TRUE)) %>%
  full_join(plotDat_long_2seasons %>% group_by(Assay) %>%
              summarise(negative = sum(statistic < 0, na.rm = TRUE))) %>%
  mutate(select = ifelse(postitive > 8*0.75 | negative > 8*0.75, TRUE, NA)) %>%
  filter(select == TRUE)

picked_proteins <- ZirFlu$protein_annot %>% filter(Assay %in% highlight_proteins$Assay)
selected_proteins <- picked_proteins$OlinkID

# metabolite with lm() model ---------------------------------------------------
inputDat_metabolite <- ZirFlu$donorSamples %>% 
  full_join(ZirFlu$donorInfo) %>% filter(time == "T1") %>%
  filter(patientID %in% ZirFlu$HAItiter$patientID) %>%
  left_join(ZirFlu$metabolite_dat %>% rownames_to_column("probenID"))

# per season
input_metabolite <- list()
input_metabolite[["2019"]] <- ZirFlu$HAItiter_2019 %>% 
  left_join(inputDat_metabolite %>% filter(season == "2019"))
input_metabolite[["2020"]] <- ZirFlu$HAItiter_2020 %>% 
  left_join(inputDat_metabolite %>% filter(season == "2020"))

lmRes_metabolite <- list()
for (season in c("2019", "2020")) {
  lmRes_metabolite[[season]]$lmRes <- get.lm_forDat(input_metabolite[[season]], 
                                                    varList = names(ZirFlu$metabolite_dat))
  lmRes_metabolite[[season]]$lmAdj <- lmRes_metabolite[[season]]$lmRes %>% get.lm_padj()
  lmRes_metabolite[[season]]$sig <- lmRes_metabolite[[season]]$lmRes %>% get.var_associFactor()
  lmRes_metabolite[[season]]$sigAdj <- lmRes_metabolite[[season]]$lmAdj %>% get.var_associFactor_padj()
}

# save result
save(lmRes_metabolite, file = "temp/lmResult_metabolite_2seasons.RData")
## venn diagram (need to check later) -------------------
load("temp/lmResult_metabolite_2seasons.RData")

## check top metabolite ----------------------------
metaSig_onlyAb_2019 <- unique(unlist(lmRes_metabolite$`2019`$sig$var_onlyAbTiter))
metaSig_onlyAb_2020 <- unique(unlist(lmRes_metabolite$`2020`$sig$var_onlyAbTiter))

metaSig_onlyAb_padj_2019 <- unique(unlist(lmRes_metabolite$`2019`$sigAdj$var_onlyAbTiter))
metaSig_onlyAb_padj_2020 <- unique(unlist(lmRes_metabolite$`2020`$sigAdj$var_onlyAbTiter))

a <- get.metaboliteAnnot(ZirFlu$metabolite_annot, c(metaSig_onlyAb_padj_2019, metaSig_onlyAb_padj_2020)) 
# heatmap --------------------
selected_metabolites <- c(metaSig_onlyAb_padj_2019, metaSig_onlyAb_padj_2020)
selected_metabolites <- c(metaSig_onlyAb_2019, metaSig_onlyAb_2020)

plotDat_wide <- list()
plotDat_long <- list()
for (season in c("2019", "2020")) {
  varlmAbTiter <- lmRes_metabolite[[season]]$lmRes %>% 
    lapply(function(x) x%>% filter(independentVariable == "abTiter")) %>%
    get.lmStatistic()
  
  heatmapDat <- varlmAbTiter[which(rownames(varlmAbTiter) %in% selected_metabolites),]
  
  plotDat_wide[[season]] <- heatmapDat %>% #as.data.frame %>% 
    rownames_to_column("ionIdx") %>% mutate(ionIdx = as.numeric(ionIdx)) %>%
    left_join(ZirFlu$metabolite_annot %>% select(ionIdx, Formula) %>% distinct()) %>%
    select(-c(ionIdx)) %>% column_to_rownames("Formula")
  
  plotDat_long[[season]] <- heatmapDat %>% #as.data.frame() %>% 
    rownames_to_column("ionIdx") %>% mutate(ionIdx = as.numeric(ionIdx)) %>%
    left_join(ZirFlu$metabolite_annot %>% select(ionIdx, Formula) %>% distinct()) %>%
    pivot_longer(cols = c(2:5), names_to = "time", values_to = "statistic") %>%
    get.dat_sigPvalue_metabolite(., lm_sigPvalue = lmRes_metabolite[[season]]$sig$var_onlyAbTiter) %>%
    rename("significance_pvalue" = "significance") %>%
    get.dat_sigPvalue_metabolite(., lm_sigPvalue = lmRes_metabolite[[season]]$sigAdj$var_onlyAbTiter)
}

plotDat_wide_2seasons <- plotDat_wide$`2019` %>% rownames_to_column("metabolite") %>%
  rename("2019_H1N1" = "H1N1_ab", "2019_H3N2" = "H3N2_ab",
         "2019_Bvic" = "Bvictor", "2019_Byam"= "Byamaga") %>% 
  full_join(plotDat_wide$`2020` %>% rownames_to_column("metabolite") %>%
              rename("2020_H1N1" = "H1N1_ab", "2020_H3N2" = "H3N2_ab",
                     "2020_Bvic" = "Bvictor", "2020_Byam" = "Byamaga")) %>%
  column_to_rownames("metabolite")

plotDat_wide_2seasons %>% as.matrix() %>% heatmap(Colv = NA, Rowv = NA)
plotDat_wide_2seasons %>% as.matrix() %>% Heatmap()
plotDat_wide_2seasons %>% as.matrix() %>% Heatmap(cluster_columns = FALSE)

plotDat_long_2seasons <- plotDat_long$`2019` %>% mutate(season = "2019") %>%
  rbind(plotDat_long$`2020` %>% mutate(season = "2020")) %>%
  unite(time, season, time) %>% mutate(time = substring(time, 1, 9))

# # if need heat map with consistent within strain  across season
# plotDat_long_2seasons$time <- factor(plotDat_long_2seasons$time,
#                                      c("2019_H1N1", "2020_H1N1",
#                                        "2019_H3N2", "2020_H3N2",
#                                        "2019_Bvic", "2020_Bvic",

plotDat_long_2seasons  %>% 
  ggplot(aes(x = time, y = Formula)) + 
  geom_tile(aes(fill = statistic)) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plotDat_long_2seasons %>% 
  ggplot(aes(x = time, y = Formula)) + 
  geom_tile(aes(fill = statistic)) +
  scale_y_discrete(limits = plotDat_long_2seasons$Formula[hclust(dist(plotDat_wide_2seasons))$order])+
  geom_text(aes(label = ifelse(significance == TRUE, "**", ""))) +
  geom_text(aes(label = ifelse(significance_pvalue == TRUE, "*", ""))) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") + theme_bw() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

# only select metabolites show consistent trend (>75% across strains and years)
highlight_metabolites <- plotDat_long_2seasons %>% group_by(Formula) %>%
  summarise(postitive = sum(statistic > 0, na.rm = TRUE)) %>%
  full_join(plotDat_long_2seasons %>% group_by(Formula) %>%
              summarise(negative = sum(statistic < 0, na.rm = TRUE))) %>%
  mutate(select = ifelse(postitive >= 8*0.75 | negative >= 8*0.75, TRUE, NA)) %>%
  filter(select == TRUE)

picked_metabolites <- ZirFlu$metabolite_annot %>% filter(Formula %in% highlight_metabolites$Formula)
selected_metabolites <- unique(picked_metabolites$ionIdx)

load("temp/20221222_keggID_allMetabolites.RData")
selected_metabolites_T1 <- keggID_allMetabolites %>% filter(ionIdx %in% selected_metabolites)

write.table(selected_metabolites_T1, file = "temp/20221222_selected_metabolites_T1.txt", sep = "\t",
            col.names = TRUE, row.names = FALSE, quote = FALSE)
write.table(keggID_allMetabolites$cpdId, file = "temp/20221222_keggID_allMetabolites.txt", sep = "\t",
            col.names = TRUE, row.names = FALSE, quote = FALSE)
