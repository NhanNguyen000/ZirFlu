library(fgsea)
library(tidyverse)

load("temp/lmResult_metabolite_2seasons.RData")
load("temp/20221019_keggIDusingFormula_library.RData")

# check the possible pathway -----------------------
sigMeta_inPathway <- list()
for (season in names(lmRes_metabolite)) {
  for (type in names(lmRes_metabolite[[season]]$sig$var_onlyAbTiter)) {
    pathway_temp <- c()
    for (i in 1:nrow(keggID_libraryIonId)) {
      meta_selected <- lmRes_metabolite[[season]]$sig$var_onlyAbTiter[[type]]
      pathways <- keggID_libraryIonId$ionIdx[i] %>% strsplit(split = "; ") %>% unlist()
      if (length(intersect(meta_selected, pathways)) > 0) {
        dat_temp <- c(keggID_libraryIonId$pathway[i],
                      length(intersect(meta_selected, pathways))) %>%
          as.data.frame() %>% t()
        pathway_temp <- rbind(pathway_temp, dat_temp)
      }
    }
    if (length(pathway_temp) > 0) {
      sigMeta_inPathway[[season]][[type]] <- pathway_temp %>%
        as.data.frame() %>% rename("pathway" = 1, "count" = 2)
    }
  }
}

sigMeta_inPathway_2019 <- sigMeta_inPathway$`2019` %>% 
  lapply(function(x) x %>% select(pathway)) %>% reduce(inner_join)

sigMeta_inPathway_2020 <- sigMeta_inPathway$`2020` %>% 
  lapply(function(x) x %>% select(pathway)) %>% reduce(inner_join)

a <- intersect(sigMeta_inPathway_2019, sigMeta_inPathway_2020) # 7 metabolite pathways
b <- keggID_libraryIonId %>% filter(pathway %in% a$pathway)

# pahtway analysis with fgsea  -----------------------
# NOTE: metabolite set is duplicate, need to re run
metabolite_sets <- list()
for(i in 1:nrow(keggID_libraryIonId)) {
  metabolite_sets[[keggID_libraryIonId$pathway[i]]] <- keggID_libraryIonId$ionIdx[i] %>%
    strsplit(split = "; ") %>% unlist() %>% unique()
}


metabolite_bigSets <- list()
for (pathway in names(metabolite_sets)) {
  if (length(metabolite_sets[[pathway]]) >=2) {
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

# focus on pathway analysis with sig. metabolite (need to check latter)---------------------
gsea_pathways[1:4]
selected_gsea_pathways <- list()
for(type in names(gsea_pathways)) {
  idx <- which(gsea_pathways[[type]]$pathway %in% sigMeta_inPathway[[type]]$pathway)
  selected_gsea_pathways[[type]] <- gsea_pathways[[type]] %>%
    slice(idx)
}

selected_gsea_pathways$T1_H1N1_titerFC %>% 
  mutate(pathway = substring(pathway, 1, 50)) %>%
  ggplot(aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=pval<0.05)) + coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA")

selected_pathway_count <- list()
for (type in names(sigMeta_inPathway)) {
  selected_pathway_count[[type]] <- sigMeta_inPathway[[type]] %>% 
    filter(count > 1)
}

selected_gsea_pathways2 <- list()
for(type in names(gsea_pathways)) {
  idx <- which(gsea_pathways[[type]]$pathway %in% sigMeta_inPathway[[type]]$pathway)
  selected_gsea_pathways2[[type]] <- gsea_pathways[[type]] %>%
    slice(idx)
}

selected_gsea_pathways2$T1_H1N1_titerFC %>% 
  mutate(pathway = substring(pathway, 1, 50)) %>%
  ggplot(aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=pval<0.05)) + coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA")

selected_gsea_pathways2$T1_Byamagata.Phuket_titerFC %>% 
  mutate(pathway = substring(pathway, 1, 50)) %>%
  ggplot(aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=pval<0.05)) + coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA")

conditions <- names(gsea_pathways)[1:4]
plotList <- list()
for (type in conditions) {
  plotDat <- selected_gsea_pathways2[[type]] %>% 
    mutate(pathway = substring(pathway, 1, 50))
  plotList[[type]] <- ggplot(plotDat, aes(reorder(pathway, NES), NES)) +
    geom_col(aes(fill=padj<0.05)) + coord_flip() +
    labs(x="Pathway", y="Normalized Enrichment Score",
         title="Hallmark pathways NES from GSEA") + 
    theme_minimal()
}
cowplot::plot_grid(plotList$H, plotList$T1_H3N2_titerFC,
                   plotList$T1_Bvictoria.Maryland_titerFC,
                   plotList$T1_Byamagata.Phuket_titerFC)

intersect(pathway$T1_H1N1_titerFC, pathway$T1_H3N2_titerFC)
intersect(pathway$T1_Bvictoria.Maryland_titerFC, pathway$T1_Byamagata.Phuket_titerFC)

# apply the heatmap -------------------------------------------------------
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

NES <- list()
for (season in names(gsea_pathways_bigSet)) {
  NES[[season]] <- gsea_pathways_bigSet[[season]] %>% 
    lapply(function(x) x %>% select(pathway, NES)) %>% 
    bind_rows(.id = "groups") %>% mutate(groups = substring(groups, 1, 7)) %>%
    pivot_wider(names_from = groups, values_from = NES)
}

plotDat_wide <- NES$`2019` %>% rename_at(c(2:5), ~paste0("2019_", .x)) %>%
  full_join(NES$`2020` %>% rename_at(c(2:5), ~paste0("2020_", .x))) %>%
  #slice(which(pathway %in% c(sigAdj_pathway_2019, sigAdj_pathway_2020))) %>%
  filter(pathway %in% c(sig_pathway_2019, sig_pathway_2020)) %>%
  rename_at(c(2:9), ~substr(.x, 1, 9)) %>% 
  column_to_rownames("pathway") %>% as.matrix()

library(ComplexHeatmap)
plotDat_wide %>% heatmap(Colv = NA, Rowv = NA)
plotDat_wide %>% as.matrix() %>% Heatmap()
plotDat_wide %>% as.matrix() %>% Heatmap(cluster_columns = FALSE)

plotDat_long <- NES$`2019` %>%
  pivot_longer(cols = c(2:5), names_to = "time", values_to = "NES") %>%
  #get.dat_sigPvalue_pathway(., sigPval_pathway = pathway_padj_bigSet$`2019`) %>%
  get.dat_sigPvalue_pathway(., sigPval_pathway = pathway_pval_bigSet$`2019`) %>%
  mutate(time = paste0("2019_", time)) %>%
  full_join(NES$`2020` %>%
              pivot_longer(cols = c(2:5), names_to = "time", values_to = "NES") %>%
              #get.dat_sigPvalue_pathway(., sigPval_pathway = pathway_padj_bigSet$`2020`) %>%
              get.dat_sigPvalue_pathway(., sigPval_pathway = pathway_pval_bigSet$`2020`) %>%
              mutate(time = paste0("2020_", time))) %>%
  #filter(pathway %in% c(sigAdj_pathway_2019, sigAdj_pathway_2020)) %>% 
  filter(pathway %in% c(sig_pathway_2019, sig_pathway_2020)) %>%
  mutate(pathway = substring(pathway, 1, 50), time = substring(time, 1, 9))

plotDat_long  %>% 
  ggplot(aes(x = time, y = pathway)) + 
  geom_tile(aes(fill = NES)) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plotDat_long%>% 
  ggplot(aes(x = time, y = pathway)) + 
  geom_tile(aes(fill = NES)) +
  scale_y_discrete(limits = plotDat_long$pathway[hclust(dist(plotDat_wide))$order])+
  geom_text(aes(label = ifelse(significance == TRUE, "*", ""))) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme_bw()

# interested pathway (need to look at later) --------------------------
upReg_pathways <- heatmapDat$pathway[c(9, 4, 25, 39, 16, 40, 27, 13, 11, 47, 41, 18)]
dowReg_pathways <- heatmapDat$pathway[c(28, 15, 43, 30, 12, 14)]

metaPath_select <- metabolite_sets[c(upReg_pathways, dowReg_pathways)]
a <- sigMeta_inPathway[1:4] %>% 
  lapply(function(x) x %>% select(pathway)) %>%
  reduce(full_join)
a <- sigMeta_inPathway[1:4] %>% 
  reduce(full_join, by = c("pathway"))
names(a)[2:5] <- names( sigMeta_inPathway)[1:4]
b <- a %>% filter(pathway %in% names(metaPath_select))

intersect(
  metaPath_select$`C-type lectin receptor signaling pathway`,
  metabolite_sig$var_onlyAbTiter$T1_H1N1_titerFC
)
intersect(
  metaPath_select$`Caffeine metabolism`,
  metabolite_sig$var_onlyAbTiter$T1_H1N1_titerFC
)

intersect(
  metaPath_select$`Biosynthesis of various alkaloids; Including: Cucurbitacin biosynthesis, Solanine and tomatine biosynthesis, Ephedrine biosynthesis, Capsaicin biosynthesis, Acridone alkaloid biosynthesis`,
  metabolite_sig$var_onlyAbTiter$T1_H1N1_titerFC
)

intersect(
  metaPath_select$`Cutin, suberine and wax biosynthesis`,
  metabolite_sig$var_onlyAbTiter$T1_Bvictoria.Maryland_titerFC
)

intersect(
  metaPath_select$`Hepatocellular carcinoma`,
  metabolite_sig$var_onlyAbTiter$T1_Bvictoria.Maryland_titerFC
)

intersect(
  metaPath_select$`Butanoate metabolism`,
  metabolite_sig$var_onlyAbTiter$T1_H1N1_titerFC
)

intersect(
  metaPath_select$`Vitamin digestion and absorption`,
  metabolite_sig$var_onlyAbTiter$T1_Byamagata.Phuket_titerFC
)

intersect(
  metaPath_select$`Sesquiterpenoid and triterpenoid biosynthesis`,
  metabolite_sig$var_onlyAbTiter$T1_H3N2_titerFC
)
# check if the pathway have only 1 metabolite --> load

