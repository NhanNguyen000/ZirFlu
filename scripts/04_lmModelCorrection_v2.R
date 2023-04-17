library(tidyverse)
library(lme4)
library(lmerTest)
library(broom.mixed)
library(ggpubr)
library(ggVennDiagram)
library(ComplexHeatmap)
# specific function for this code ---------------------------
get.lm_forDat <- function(inputDat, varList) {
  lm_results <- list()
  timeList <- c("T1", "T3", "T4")
  abTiters <- c("H1N1_titerFC", "H3N2_titerFC", 
                "Bvictoria.Maryland_titerFC", "Byamagata.Phuket_titerFC")
  for (timepoint in timeList) {
    dat_time <- inputDat %>% filter(time %in% timepoint)
    for (name in abTiters) {
      resName <- paste0(timepoint, "_", name)
      lm_results[[resName]] <- get.lmTest_correction(variableList = varList,
                                                     abTiter = name, inputDat = dat_time)
    }
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
get.var_associFactor_padj <- function(lm_results) {
  var_associDisease <- list()
  var_associAbTiter <- list()
  var_onlyAbTiter <- list()
  
  for (i in names(lm_results)) {
    var_associDisease[[i]] <- lm_results[[i]] %>% filter(p.adj < 0.05) %>%
      filter(independentVariable == "diseasehealthy")
    
    var_associAbTiter[[i]] <- lm_results[[i]] %>% filter(p.adj < 0.05) %>%
      filter(independentVariable == "abTiter")
    
    var_onlyAbTiter[[i]] <- get.var_onlyAbTiter(var_associAbTiter, 
                                                var_associDisease,
                                                abTiter = i)
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
# venn diagram
get.vennDat <- function(inputDat, time) {
  dat_temp <- inputDat %>% lapply(function(x) x%>% select(targetVariable))
  outcome <- dat_temp[grep(time, names(dat_temp))] %>% unlist(recursive = FALSE)
  names(outcome) <- substring(names(outcome), 1, 7)
  return(outcome)
}

get.vennPlot_pertime <- function(var_associDisease, var_associAbTiter, time) {
  cowplot::plot_grid(
    ggVennDiagram(get.vennDat(var_associDisease, time),
                  label_alpha = 0) + scale_fill_gradient(low="white",high = "blue"),
    ggVennDiagram(get.vennDat(var_associAbTiter, time),
                  label_alpha = 0) + scale_fill_gradient(low="white",high = "blue"),
    labels = c("Disease-associated", "abTiterFC-associated"), nrow = 2)
}

get.vennPlot_allTime <- function(var_associDisease, var_associAbTiter) {
  cowplot::plot_grid(
    ggVennDiagram(get.vennDat(var_associDisease, "T1"),
                  label_alpha = 0) + scale_fill_gradient(low="white",high = "blue"),
    ggVennDiagram(get.vennDat(var_associAbTiter, "T1"),
                  label_alpha = 0) + scale_fill_gradient(low="white",high = "blue"),
    ggVennDiagram(get.vennDat(var_associDisease, "T3"),
                  label_alpha = 0) + scale_fill_gradient(low="white",high = "blue"),
    ggVennDiagram(get.vennDat(var_associAbTiter, "T3"),
                  label_alpha = 0) + scale_fill_gradient(low="white",high = "blue"),
    ggVennDiagram(get.vennDat(var_associDisease, "T4"),
                  label_alpha = 0) + scale_fill_gradient(low="white",high = "blue"),
    ggVennDiagram(get.vennDat(var_associAbTiter, "T4"),
                  label_alpha = 0) + scale_fill_gradient(low="white",high = "blue"),
    byrow = FALSE
  )
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
      mutate(significance = ifelse(OlinkID %in% lm_sigPvalue[[name]], TRUE, NA))
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
# protein -----------------
inputDat_protein <- ZirFlu$donorInfo %>% full_join(ZirFlu$donorSamples) %>%
  select(patientID, probenID, sex, age, condition, category, time) %>%
  full_join(ZirFlu$HAItiter_Feb2022 %>% select(-condition)) %>%
  right_join(ZirFlu$protein_dat %>% rownames_to_column("probenID") %>% 
               mutate(probenID = as.numeric(probenID)))

lm_results_protein <- get.lm_forDat(inputDat_protein, varList = names(ZirFlu$protein_dat))
lm_adj_protein <- get.lm_padj(lm_results_protein)
proteins_sig <- get.var_associFactor(lm_results_protein)
proteins_sig_adj <- get.var_associFactor_padj(lm_adj_protein)

save(lm_results_protein, lm_adj_protein, proteins_sig, proteins_sig_adj,
     file = "temp/lmResult_protein.RData")

## venn diagram -------------------
load("temp/lmResult_protein.RData")
# at T1
get.vennPlot_pertime(var_associDisease = proteins_sig$var_associDisease, 
                     var_associAbTiter = proteins_sig$var_associAbTiter,
                     time = "T1")
get.vennPlot_pertime(var_associDisease = proteins_sig_adj$var_associDisease, 
                     var_associAbTiter = proteins_sig_adj$var_associAbTiter,
                     time = "T1")
# in all time point
get.vennPlot_allTime(var_associDisease = proteins_sig$var_associDisease, 
                     var_associAbTiter = proteins_sig$var_associAbTiter)
get.vennPlot_allTime(var_associDisease = proteins_sig_adj$var_associDisease, 
                     var_associAbTiter = proteins_sig_adj$var_associAbTiter)

## check top protein ---
unique(unlist(proteins_sig_adj$var_onlyAbTiter))
a <- get.proteinAnnot(ZirFlu$protein_annot, 
                      unique(unlist(proteins_sig_adj$var_onlyAbTiter))) 

proteins_sigAssociDisease <- get.targetedVar(proteins_sig$var_associDisease)
length(unique(proteins_sigAssociDisease[1:4] %>% unlist()))

proteins_sigAssociAbFC <- get.targetedVar(proteins_sig$var_associAbTiter)
sigProteinT1_Ab_pvalue <- unique(proteins_sigAssociAbFC[1:4] %>% unlist())

proteins_sigAssociAbFC_padj <- get.targetedVar(proteins_sig_adj$var_associAbTiter)
sigProteinT1_Ab_padj <- unique(proteins_sigAssociAbFC_padj[1:4] %>% unlist())

sigProteinT1_onlyAb_pvalue <- unique(proteins_sig$var_onlyAbTiter[1:4] %>% unlist())
sigProteinT1_onlyAb_padj <- unique(proteins_sig_adj$var_onlyAbTiter[1:4] %>% unlist())

# heatmap --------------------
varlmAbTiter <- lm_results_protein %>% 
  lapply(function(x) x%>% filter(independentVariable == "abTiter")) %>%
  get.lmStatistic() %>% as.matrix()

# heatmapDat <- varlmAbTiter[which(rownames(varlmAbTiter) %in% unique(unlist(proteins_sig$var_onlyAbTiter))),]
heatmapDat <- varlmAbTiter[which(rownames(varlmAbTiter) %in% unique(unlist(proteins_sig_adj$var_onlyAbTiter))),]
heatmapDat <- varlmAbTiter[which(rownames(varlmAbTiter) %in% unique(unlist(proteins_sig_adj$var_onlyAbTiter[1:4]))),]
heatmapDat <- varlmAbTiter[which(rownames(varlmAbTiter) %in% unique(unlist(proteins_sig_adj$var_onlyAbTiter[1:4]))),] %>%
  as.data.frame %>% select(starts_with("T1"))

plotDat_wide <- heatmapDat %>% as.data.frame %>% 
  rownames_to_column("OlinkID") %>% left_join(ZirFlu$protein_annot) %>%
  select(-c(OlinkID, UniProt)) %>% column_to_rownames("Assay") %>% as.matrix

plotDat_wide %>% heatmap(Colv = NA, Rowv = NA)
plotDat_wide %>% as.matrix() %>% Heatmap()
plotDat_wide %>% as.matrix() %>% Heatmap(cluster_columns = FALSE)

plotDat_long <- heatmapDat %>% as.data.frame() %>% 
  rownames_to_column("OlinkID") %>% left_join(ZirFlu$protein_annot) %>%
  pivot_longer(cols = starts_with("T"), 
               names_to = "time", values_to = "statistic") %>%
  #get.dat_sigPvalue(., lm_sigPvalue = proteins_sig$var_onlyAbTiter)
  get.dat_sigPvalue(., lm_sigPvalue = proteins_sig_adj$var_onlyAbTiter)

plotDat_long  %>% 
  ggplot(aes(x = time, y = Assay)) + 
  geom_tile(aes(fill = statistic)) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plotDat_long %>% 
  ggplot(aes(x = time, y = Assay)) + 
  geom_tile(aes(fill = statistic)) +
  scale_y_discrete(limits = plotDat_long$Assay[hclust(dist(plotDat_wide))$order])+
  geom_text(aes(label = ifelse(significance == TRUE, "*", ""))) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme_bw()

# check complex heatmap package

# metabolite -----------------
inputDat_metabolite <- ZirFlu$donorInfo %>% full_join(ZirFlu$donorSamples) %>%
  select(patientID, probenID, sex, age, condition, category, time) %>%
  full_join(ZirFlu$HAItiter_Feb2022 %>% select(-condition)) %>%
  right_join(ZirFlu$metabolite_dat %>% rownames_to_column("probenID") %>% 
               mutate(probenID = as.numeric(probenID)))

lm_results_metabolite <- get.lm_forDat(inputDat_metabolite, varList = names(ZirFlu$metabolite_dat))
lm_adj_metabolite <- get.lm_padj(lm_results_metabolite)
metabolite_sig <- get.var_associFactor(lm_results_metabolite)
metabolite_sig_adj <- get.var_associFactor_padj(lm_adj_metabolite)

save(lm_results_metabolite, lm_adj_metabolite, metabolite_sig, metabolite_sig_adj,
     file = "temp/lmResult_metabolites.RData")

## venn diagram -------------------
load("temp/lmResult_metabolites.RData")

get.vennPlot_pertime(var_associDisease = metabolite_sig$var_associDisease, 
                     var_associAbTiter = metabolite_sig$var_associAbTiter,
                     time = "T1")
get.vennPlot_pertime(var_associDisease = metabolite_sig_adj$var_associDisease, 
                     var_associAbTiter = metabolite_sig_adj$var_associAbTiter,
                     time = "T1")

get.vennPlot_allTime(var_associDisease = metabolite_sig$var_associDisease, 
                     var_associAbTiter = metabolite_sig$var_associAbTiter)
get.vennPlot_allTime(var_associDisease = metabolite_sig_adj$var_associDisease, 
                     var_associAbTiter = metabolite_sig_adj$var_associAbTiter)

## check top metabolite ----------------------------
unique(unlist(metabolite_sig_adj$var_onlyAbTiter))
a <- get.metaboliteAnnot(ZirFlu$metabolite_annot, 
                      unique(unlist(metabolite_sig_adj$var_onlyAbTiter))) 

metabolite_sigAssociDisease <- get.targetedVar(metabolite_sig$var_associDisease)
metaboliteT1_disease_pvalue <- unique(metabolite_sigAssociDisease[1:4] %>% unlist())

metabolite_sigAssociAbFC <- get.targetedVar(metabolite_sig$var_associAbTiter)
metaboliteT1_Ab_pvalue <- unique(metabolite_sigAssociAbFC[1:4] %>% unlist())

metabolite_sigAssociDisease_adj <- get.targetedVar(metabolite_sig_adj$var_associDisease)
metaboliteT1_disease_padj <- unique(metabolite_sigAssociDisease_adj[1:4] %>% unlist())

metabolite_sigAssociAbFC_padj <- get.targetedVar(metabolite_sig_adj$var_associAbTiter)
metaboliteT1_Ab_padj <- unique(metabolite_sigAssociAbFC_padj[1:4] %>% unlist())

metaboliteT1_onlyAb_pvalue <- unique(metabolite_sig$var_onlyAbTiter[1:4] %>% unlist())
metaboliteT1_onlyAb_padj <- unique(metabolite_sig_adj$var_onlyAbTiter[1:4] %>% unlist())

# heatmap --------------------
varlmAbTiter <- lm_results_metabolite %>% 
  lapply(function(x) x%>% filter(independentVariable == "abTiter")) %>%
  get.lmStatistic() %>% as.matrix()

heatmapDat <- varlmAbTiter[which(rownames(varlmAbTiter) %in% unique(unlist(metabolite_sig$var_onlyAbTiter))),]
# heatmapDat <- varlmAbTiter[which(rownames(varlmAbTiter) %in% unique(unlist(metabolite_sig_adj$var_onlyAbTiter))),]
heatmapDat <- varlmAbTiter[which(rownames(varlmAbTiter) %in% unique(unlist(metabolite_sig_adj$var_onlyAbTiter[1:4]))),] %>%
  as.data.frame %>% select(starts_with("T1"))

plotDat_wide <- heatmapDat %>% as.data.frame %>% 
  rownames_to_column("ionIdx") %>% mutate(ionIdx = as.numeric(ionIdx)) %>%
  left_join(ZirFlu$metabolite_annot %>% select(ionIdx, Formula) %>% distinct()) %>%
  select(-c(ionIdx)) %>% column_to_rownames("Formula") %>% as.matrix

plotDat_wide %>% heatmap(Colv = NA, Rowv = NA)
plotDat_wide %>% as.matrix() %>% Heatmap()
plotDat_wide %>% as.matrix() %>% Heatmap(cluster_columns = FALSE)

plotDat_long <- heatmapDat %>% as.data.frame() %>% 
  rownames_to_column("ionIdx") %>% mutate(ionIdx = as.numeric(ionIdx)) %>%
  left_join(ZirFlu$metabolite_annot %>% select(ionIdx, Formula) %>% distinct()) %>%
  pivot_longer(cols = starts_with("T"), 
               names_to = "time", values_to = "statistic") %>%
  #get.dat_sigPvalue_metabolite(., lm_sigPvalue = metabolite_sig$var_onlyAbTiter)
  get.dat_sigPvalue_metabolite(., lm_sigPvalue = metabolite_sig_adj$var_onlyAbTiter)

plotDat_long  %>% 
  ggplot(aes(x = time, y = Formula)) + 
  geom_tile(aes(fill = statistic)) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plotDat_long %>% 
  ggplot(aes(x = time, y = Formula)) + 
  geom_tile(aes(fill = statistic)) +
  scale_y_discrete(limits = plotDat_long$Formula[hclust(dist(plotDat_wide))$order])+
  geom_text(aes(label = ifelse(significance == TRUE, "*", ""))) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme_bw()

# mlic and citric acid expression
citric - C6H8O7
malate - C4H6O5 - ionId 132

## T1 metabolite plot -------------------
selectedPro <- "132"
abTiter <- "H1N1_titerFC"

plot_dat <- inputDat %>% filter(time == "T1") %>% 
  select(patientID, condition, disease, category, c( "H1N1_titerFC", "132")) %>%
  rename(ionId132 = "132")

ggplot(plot_dat, aes(x = H1N1_titerFC, y = ionId132)) +
  geom_point(aes(color = category, shape = condition), 
             size = 3, alpha = 0.8,
             position = position_jitter(w = 0.08, h = 0.1)) +
  geom_smooth(method = "lm", se = FALSE) + theme_bw()

plot_dat$category <- factor(plot_dat$category, levels = c("NR", "Other", "TR"))
plot_dat$disease <- factor(plot_dat$disease, 
                           levels = c("healthy", "cirrhosis"))
plot_dat$condition <- factor(plot_dat$condition, 
                             levels = c("healthy", "compensated cirrhosis", "decompensated cirrhosis"))
cowplot::plot_grid(
  ggboxplot(plot_dat, x = "category", y = "ionId132", palette = "jco", add = "jitter"),
  ggboxplot(plot_dat, x = "disease", y = "ionId132", palette = "jco", add = "jitter"),
  ggboxplot(plot_dat, x = "condition", y = "ionId132", palette = "jco", add = "jitter") + 
    rotate_x_text(20), nrow = 1)
