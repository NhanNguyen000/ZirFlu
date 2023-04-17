library(fgsea)
library(tidyverse)

load("temp/lmResult_metabolites.RData")
load("temp/20221019_keggIDusingFormula_library.RData")

# check the possible pathway -----------------------
sigMeta_inPathway <- list()
for (type in names(metabolite_sig$var_onlyAbTiter)) {
  pathway_temp <- c()
  for (i in 1:nrow(keggID_libraryIonId)) {
    meta_selected <- metabolite_sig$var_onlyAbTiter[[type]]
    meta_pathways <- keggID_libraryIonId$ionIdx[i] %>% 
      strsplit(split = "; ") %>% unlist()
    
    if (length(intersect(meta_selected, meta_pathways)) >0) {
      dat_temp <- c(keggID_libraryIonId$pathway[i], length(intersect(meta_selected, meta_pathways))) %>%
        as.data.frame() %>% t()
      pathway_temp <- rbind(pathway_temp, dat_temp )
    }
  }
  sigMeta_inPathway[[type]] <- pathway_temp %>% 
    as.data.frame() %>% rename("pathway" = 1, "count" = 2)
}

sigMeta_inPathway_T1 <- sigMeta_inPathway[1:4] %>% 
  lapply(function(x) x %>% select(pathway)) %>%
  reduce(inner_join)


# pahtway analysis with fgsea  -----------------------
# NOTE: metabolite set is duplicate, need to re run
metabolite_sets <- list()
for(i in 1:nrow(keggID_libraryIonId)) {
  metabolite_sets[[keggID_libraryIonId$pathway[i]]] <- keggID_libraryIonId$ionIdx[i] %>%
    strsplit(split = "; ") %>% unlist() %>% unique()
}

gsea_pathways <- list()
for(type in names(lm_results_metabolite)) {
  inputDat_temp <- lm_results_metabolite[[type]] %>% 
    filter(independentVariable == "abTiter")
  
  ranked_list <- inputDat_temp$statistic
  names(ranked_list) <- inputDat_temp$targetVariable
  sort(ranked_list, decreasing = T) -> ranked_list
  
  gsea_pathways[[type]] <- fgsea(pathways = metabolite_sets, 
                    stats = ranked_list, nPermSimple = 5000)
}

plotList <- list()
for (type in names(gsea_pathways)) {
  plotDat <- gsea_pathways[[type]] %>% filter(pval <0.05) %>% 
    mutate(pathway = substring(pathway, 1, 50))
  plotList[[type]] <- ggplot(plotDat, aes(reorder(pathway, NES), NES)) +
    geom_col(aes(fill=padj<0.05)) + coord_flip() +
    labs(x="Pathway", y="Normalized Enrichment Score",
         title="Hallmark pathways NES from GSEA") + 
    theme_minimal()
}
cowplot::plot_grid(plotList$T1_H1N1_titerFC, plotList$T1_H3N2_titerFC,
                   plotList$T1_Bvictoria.Maryland_titerFC,
                   plotList$T1_Byamagata.Phuket_titerFC)

# focus on T1 pathway analysis ---------------------
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
cowplot::plot_grid(plotList$T1_H1N1_titerFC, plotList$T1_H3N2_titerFC,
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

get.sig_padj <- function(lm_results) {
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

pathway_pval <- get.sig_pval(gsea_pathways)
#pathway_padj <- get.sig_padj(gsea_pathways)

sigPathway_T1 <- unique(pathway_pval[1:4] %>% 
                          lapply(function(x) x %>% select(pathway)) %>% unlist())

NES <- gsea_pathways %>% lapply(function(x) x %>% select(pathway, NES)) %>% 
  bind_rows(.id = "groups") %>% mutate(groups = substring(groups, 1, 7)) %>%
  pivot_wider(names_from = groups, values_from = NES)
  

heatmapDat <- NES[which(NES$pathway %in% sigPathway_T1),] %>%
  as.data.frame %>% select(pathway, starts_with("T1")) 

library(ComplexHeatmap)
plotDat_wide <- heatmapDat %>% column_to_rownames("pathway") %>% as.matrix()
plotDat_wide %>% heatmap(Colv = NA, Rowv = NA)
plotDat_wide %>% as.matrix() %>% Heatmap()
plotDat_wide %>% as.matrix() %>% Heatmap(cluster_columns = FALSE)

plotDat_long <- heatmapDat %>% as.data.frame() %>%
  pivot_longer(cols = starts_with("T"), 
               names_to = "time", values_to = "NES") %>%
  get.dat_sigPvalue_pathway(., sigPval_pathway = pathway_pval) %>%
  mutate(pathway = substring(pathway, 1, 50))

plotDat_long  %>% 
  ggplot(aes(x = time, y = pathway)) + 
  geom_tile(aes(fill = NES)) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plotDat_long %>% 
  ggplot(aes(x = time, y = pathway)) + 
  geom_tile(aes(fill = NES)) +
  scale_y_discrete(limits = plotDat_long$pathway[hclust(dist(plotDat_wide))$order])+
  geom_text(aes(label = ifelse(significance == TRUE, "*", ""))) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme_bw()

# interested pathway:
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

