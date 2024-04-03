library(fgsea)
library(tidyverse)
library(ComplexHeatmap)
# function ----------------------------------------------------------
get.related_pathways <- function(keggID_libraryIonId, selected_metabolites) {
  pathway_temp <- c()
  for (i in 1:nrow(keggID_libraryIonId)) {
    pathways <- keggID_libraryIonId$ionIdx[i] %>% strsplit(split = "; ") %>% unlist()
    if (length(intersect(selected_metabolites, pathways)) > 0) {
      dat_temp <- c(keggID_libraryIonId$pathway[i],
                    length(intersect(selected_metabolites, pathways))) %>%
        as.data.frame() %>% t()
      pathway_temp <- rbind(pathway_temp, dat_temp)
    }
  }
  if (length(pathway_temp) > 0) {
    outcome <- pathway_temp %>% as.data.frame() %>% rename("pathway" = 1, "count" = 2)
  } else outcome <- NULL
  return(outcome)
}

get.sig_pval <- function(inputDat) {
  sig_pval <- list()
  for (i in names(inputDat)) {
    sig_pval[[i]] <- inputDat[[i]] %>% filter(pval < 0.05)
  }
  return(sig_pval)
}

get.sig_padj <- function(inputDat) {
  sig_padj <- list()
  for (i in names(inputDat)) {
    sig_padj[[i]] <- inputDat[[i]] %>% filter(padj < 0.05)
  }
  return(sig_padj)
}

get.dat_sigPvalue_pathway <- function(inputDat, sigPval_pathway) {
  outcome <- as.data.frame(matrix(nrow=0, ncol = 4))
  for (name in names(sigPval_pathway)) {
    dat_temp <- inputDat %>% filter(time == substring(name, 1, 7)) %>% 
      mutate(significance = ifelse(pathway %in% sigPval_pathway[[name]]$pathway, TRUE, NA))
    outcome <- rbind(outcome, dat_temp)
  }
  return(outcome)
}

# T1 - check the possible pathway ----------------------------------------------------
load("temp/lmResult_metabolite_2seasons.RData")
load("temp/20221019_keggIDusingFormula_library.RData")

sigMeta_inPathway <- list()
for (season in names(lmRes_metabolite)) {
  for (type in names(lmRes_metabolite[[season]]$sig$var_onlyAbTiter)) {
    sigMeta_inPathway[[season]][[type]] <- get.related_pathways(keggID_libraryIonId, 
        selected_metabolites = lmRes_metabolite[[season]]$sig$var_onlyAbTiter[[type]])
    
  }
}

sigMeta_inPathway_2019 <- sigMeta_inPathway$`2019` %>% 
  lapply(function(x) x %>% select(pathway)) %>% reduce(inner_join)

sigMeta_inPathway_2020 <- sigMeta_inPathway$`2020` %>% 
  lapply(function(x) x %>% select(pathway)) %>% reduce(inner_join)

a <- intersect(sigMeta_inPathway_2019, sigMeta_inPathway_2020) # 7 metabolite pathways
b <- keggID_libraryIonId %>% filter(pathway %in% a$pathway)

## pathway analysis with fgsea --------------------------------------------------
metabolite_sets <- list()
for(i in 1:nrow(keggID_libraryIonId)) {
  metabolite_sets[[keggID_libraryIonId$pathway[i]]] <- keggID_libraryIonId$ionIdx[i] %>%
    strsplit(split = "; ") %>% unlist() %>% unique()
}

# add the "Prostaglandin Synthesis and Regulation" pathway based on the cpdb analysis
metabolite_sets$`Prostaglandin Synthesis and Regulation` <- c("782", "935", "944")

metabolite_bigSets <- list()
for (pathway in names(metabolite_sets)) {
  if (length(metabolite_sets[[pathway]]) >2) {
    metabolite_bigSets[[pathway]] <- metabolite_sets[[pathway]]
  }
}

gsea_pathways <- list()
gsea_pathways_bigSet <- list()
for(season in names(lmRes_metabolite)) {
  for(type in names(lmRes_metabolite[[season]]$lmRes)) {
    inputDat_temp <- lmRes_metabolite[[season]]$lmRes[[type]] %>% 
      filter(independentVariable == "abTiter")
    
    ranked_list <- inputDat_temp$statistic
    names(ranked_list) <- inputDat_temp$targetVariable
    sort(ranked_list, decreasing = T) -> ranked_list
    
    gsea_pathways[[season]][[type]] <- fgsea(pathways = metabolite_sets, 
                                             stats = ranked_list, nPermSimple = 5000)
    gsea_pathways_bigSet[[season]][[type]] <- fgsea(pathways = metabolite_bigSets, 
                                                    stats = ranked_list, nPermSimple = 5000)
  }
}


plotList <- list()
for(season in names(gsea_pathways)) {
  for (type in names(gsea_pathways[[season]])) {
    plotDat <- gsea_pathways[[season]][[type]] %>% 
      filter(pval <0.05) %>% mutate(pathway = substring(pathway, 1, 50))
    
    plotList[[season]][[type]] <- ggplot(plotDat, aes(reorder(pathway, NES), NES)) +
      geom_col(aes(fill=padj<0.05)) + coord_flip() +
      labs(x="Pathway", y="Normalized Enrichment Score",
           title="Hallmark pathways NES from GSEA") + 
      theme_minimal()
  }
}

cowplot::plot_grid(plotList$`2019`$H1N1_abFC, plotList$`2019`$H3N2_abFC,
                   plotList$`2019`$Bvictoria_abFC, plotList$`2019`$Byamagata_abFC)
## apply the heatmap -------------------------------------------------------
pathway_pval <- list()
pathway_padj <- list()
pathway_pval_bigSet <- list()
pathway_padj_bigSet <- list()
for (season in names(gsea_pathways)) {
  pathway_pval[[season]] <- get.sig_pval(gsea_pathways[[season]])
  pathway_padj[[season]] <- get.sig_padj(gsea_pathways[[season]])
  pathway_pval_bigSet[[season]] <- get.sig_pval(gsea_pathways_bigSet[[season]])
  pathway_padj_bigSet[[season]] <- get.sig_padj(gsea_pathways_bigSet[[season]])
}

sig_pathway_2019 <- unique(pathway_pval_bigSet$`2019` %>% 
                             lapply(function(x) x %>% select(pathway)) %>% unlist())
sig_pathway_2020 <- unique(pathway_pval_bigSet$`2020` %>% 
                             lapply(function(x) x %>% select(pathway)) %>% unlist())

sigAdj_pathway_2019 <- unique(pathway_padj_bigSet$`2019` %>% 
                                lapply(function(x) x %>% select(pathway)) %>% unlist())
sigAdj_pathway_2020 <- unique(pathway_padj_bigSet$`2020` %>% 
                                lapply(function(x) x %>% select(pathway)) %>% unlist())

## making heatmap ----------------------------------------------------------
selected_pathways <- c(sig_pathway_2019, sig_pathway_2020)

selected_pathways <- c("Linoleic acid metabolism",
                       "Pentose and glucuronate interconversions",
                       "Prostaglandin Synthesis and Regulation",
                       "Steroid hormone biosynthesis")

NES <- list()
plotDat_wide <- list()
plotDat_long <- list()
for (season in names(gsea_pathways_bigSet)) {
  NES[[season]] <- gsea_pathways_bigSet[[season]] %>% 
    lapply(function(x) x %>% select(pathway, NES)) %>% 
    bind_rows(.id = "groups") %>% mutate(groups = substring(groups, 1, 7)) %>%
    pivot_wider(names_from = groups, values_from = NES)
  
  plotDat_wide[[season]] <- NES[[season]] %>% filter(pathway %in% selected_pathways)
  
  plotDat_long[[season]] <- NES[[season]] %>% filter(pathway %in% selected_pathways) %>%
    pivot_longer(cols = c(2:5), names_to = "time", values_to = "NES") %>%
    get.dat_sigPvalue_pathway(., sigPval_pathway = pathway_pval_bigSet[[season]]) %>%
    rename("significance_pvalue" = "significance") %>%
    get.dat_sigPvalue_pathway(., sigPval_pathway = pathway_padj_bigSet[[season]])
}

plotDat_wide_2seasons <- plotDat_wide$`2019` %>%
  rename("2019_H1N1" = "H1N1_ab", "2019_H3N2" = "H3N2_ab",
         "2019_Bvic" = "Bvictor", "2019_Byam"= "Byamaga") %>% 
  full_join(plotDat_wide$`2020` %>%
              rename("2020_H1N1" = "H1N1_ab", "2020_H3N2" = "H3N2_ab",
                     "2020_Bvic" = "Bvictor", "2020_Byam" = "Byamaga")) %>%
  column_to_rownames("pathway") %>% as.matrix()

# plot wide data
plotDat_wide_2seasons %>% heatmap(Colv = NA, Rowv = NA)
plotDat_wide_2seasons %>% Heatmap()
plotDat_wide_2seasons %>% Heatmap(cluster_columns = FALSE)

# plot long data
plotDat_long_2seasons <- plotDat_long$`2019` %>% mutate(season = "2019") %>%
  rbind(plotDat_long$`2020` %>% mutate(season = "2020")) %>%
  unite(time, season, time) %>%
  mutate(pathway = substring(pathway, 1, 50), time = substring(time, 1, 9))

plotDat_long_2seasons %>% 
  ggplot(aes(x = time, y = pathway)) + 
  geom_tile(aes(fill = NES)) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plotDat_long_2seasons %>% 
  ggplot(aes(x = time, y = pathway)) + 
  geom_tile(aes(fill = NES)) +
  scale_y_discrete(limits = plotDat_long_2seasons$pathway[hclust(dist(plotDat_wide_2seasons))$order])+
  geom_text(aes(label = ifelse(significance == TRUE, "**", ""))) +
  geom_text(aes(label = ifelse(significance_pvalue == TRUE, "*", ""))) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") + theme_bw() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

# T2 and T3 - check the possible pathway ----------------------------------------------------
load("temp/20221019_keggIDusingFormula_library.RData")

sigMeta_inPathway <- list()
for (season in names(resPrePost_metabolite)) {
  for (type in names(resPrePost_metabolite[[season]])) {
    sigMeta_inPathway[[season]][[type]] <- get.related_pathways(keggID_libraryIonId, 
        selected_metabolites = rownames(resPrePost_metabolite[[season]][[type]]$pval))
    
  }
}

## pathway analysis with fgsea --------------------------------------------------
metabolite_sets <- list()
for(i in 1:nrow(keggID_libraryIonId)) {
  metabolite_sets[[keggID_libraryIonId$pathway[i]]] <- keggID_libraryIonId$ionIdx[i] %>%
    strsplit(split = "; ") %>% unlist() %>% unique()
}

# add the "Prostaglandin Synthesis and Regulation" pathway based on the cpdb analysis
metabolite_sets$`Prostaglandin Synthesis and Regulation` <- c("782", "935", "944")

metabolite_bigSets <- list()
for (pathway in names(metabolite_sets)) {
  if (length(metabolite_sets[[pathway]]) >2) {
    metabolite_bigSets[[pathway]] <- metabolite_sets[[pathway]]
  }
}

gsea_pathways <- list()
gsea_pathways_bigSet <- list()
for(season in names(resPrePost_metabolite)) {
  for(type in names(resPrePost_metabolite[[season]])) {
    inputDat_temp <- resPrePost_metabolite[[season]][[type]]$resTable
    
    ranked_list <- inputDat_temp$t
    names(ranked_list) <- rownames(inputDat_temp)
    sort(ranked_list, decreasing = T) -> ranked_list
    
    gsea_pathways[[season]][[type]] <- fgsea(pathways = metabolite_sets, 
                                             stats = ranked_list, nPermSimple = 5000)
    gsea_pathways_bigSet[[season]][[type]] <- fgsea(pathways = metabolite_bigSets, 
                                                    stats = ranked_list, nPermSimple = 5000)
  }
}
## apply the heatmap -------------------------------------------------------
pathway_pval <- list()
pathway_padj <- list()
pathway_pval_bigSet <- list()
pathway_padj_bigSet <- list()
for (season in names(gsea_pathways)) {
  pathway_pval[[season]] <- get.sig_pval(gsea_pathways[[season]])
  pathway_padj[[season]] <- get.sig_padj(gsea_pathways[[season]])
  pathway_pval_bigSet[[season]] <- get.sig_pval(gsea_pathways_bigSet[[season]])
  pathway_padj_bigSet[[season]] <- get.sig_padj(gsea_pathways_bigSet[[season]])
}

sig_pathway_2019 <- unique(pathway_pval_bigSet$`2019` %>% 
                             lapply(function(x) x %>% select(pathway)) %>% unlist())
sig_pathway_2020 <- unique(pathway_pval_bigSet$`2020` %>% 
                             lapply(function(x) x %>% select(pathway)) %>% unlist())

sigAdj_pathway_2019 <- unique(pathway_padj_bigSet$`2019` %>% 
                                lapply(function(x) x %>% select(pathway)) %>% unlist())
sigAdj_pathway_2020 <- unique(pathway_padj_bigSet$`2020` %>% 
                                lapply(function(x) x %>% select(pathway)) %>% unlist())

## making heatmap ----------------------------------------------------------
selected_pathways <- c(sig_pathway_2019, sig_pathway_2020)

selected_pathways <- c("Galactose metabolism",
                       "Steroid hormone biosynthesis",
                       "Fructose and mannose metabolism",
                       "Glycolysis / Gluconeogenesis")

a<- keggID_libraryIonId %>% filter(pathway %in% selected_pathways)

NES <- list()
plotDat_wide <- list()
plotDat_long <- list()
for (season in names(gsea_pathways_bigSet)) {
  NES[[season]] <- gsea_pathways_bigSet[[season]] %>% 
    lapply(function(x) x %>% select(pathway, NES)) %>% 
    bind_rows(.id = "groups") %>% mutate(groups = substring(groups, 1, 7)) %>%
    pivot_wider(names_from = groups, values_from = NES)
  
  plotDat_wide[[season]] <- NES[[season]] %>% filter(pathway %in% selected_pathways)
  
  plotDat_long[[season]] <- NES[[season]] %>% filter(pathway %in% selected_pathways) %>%
    pivot_longer(cols = c(2:3), names_to = "time", values_to = "NES") %>%
    get.dat_sigPvalue_pathway(., sigPval_pathway = pathway_pval_bigSet[[season]]) %>%
    rename("significance_pvalue" = "significance") %>%
    get.dat_sigPvalue_pathway(., sigPval_pathway = pathway_padj_bigSet[[season]])
}

plotDat_wide_2seasons <- plotDat_wide$`2019` %>%
  rename("2019_T1vsT2" = "T1vsT2", "2019_T1vsT3" = "T1vsT3") %>% 
  full_join(plotDat_wide$`2020` %>% rename("2020_T1vsT2" = "T1vsT2", "2020_T1vsT3" = "T1vsT3")) %>%
  column_to_rownames("pathway") %>% as.matrix()

# plot wide data
plotDat_wide_2seasons %>% heatmap(Colv = NA, Rowv = NA)
plotDat_wide_2seasons %>% Heatmap()
plotDat_wide_2seasons %>% Heatmap(cluster_columns = FALSE)

# plot long data
plotDat_long_2seasons <- plotDat_long$`2019` %>% mutate(season = "2019") %>%
  rbind(plotDat_long$`2020` %>% mutate(season = "2020")) %>%
  unite(time, season, time) %>%
  mutate(pathway = substring(pathway, 1, 50), time = substring(time, 1, 12))

plotDat_long_2seasons %>% 
  ggplot(aes(x = time, y = pathway)) + 
  geom_tile(aes(fill = NES)) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") + theme_bw() + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) 

plotDat_long_2seasons %>% 
  ggplot(aes(x = time, y = pathway)) + 
  geom_tile(aes(fill = NES)) +
  scale_y_discrete(limits = plotDat_long_2seasons$pathway[hclust(dist(plotDat_wide_2seasons))$order])+
  geom_text(aes(label = ifelse(significance == TRUE, "**", ""))) +
  geom_text(aes(label = ifelse(significance_pvalue == TRUE, "*", ""))) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") + theme_bw() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

# check pathways


