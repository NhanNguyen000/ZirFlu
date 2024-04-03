
get.pca_plot <- function(pca, metadat, groupType) {
  library(tidyverse)
  library(caret)
  
  ggplot(data.frame(pca$x), aes(PC1, PC2, color = metadat[[groupType]])) +
    geom_point(alpha=.5) +
    theme_classic() +
    theme(aspect.ratio = 1, plot.title = element_text(hjust = .5)) +
    ggsci::scale_color_nejm() +
    labs(x = paste0("PC1 [", 
                    round(100*summary(pca)$importance["Proportion of Variance", "PC1"], 
                          digits = 2), 
                    "%]"),
         y = paste0("PC2 [", 
                    round(100*summary(pca)$importance["Proportion of Variance", "PC2"], 
                          digits = 2), 
                    "%]"), 
         color = groupType) +
    stat_ellipse()
}

# code -------------------------------
library(tidyverse)
library(caret)
## PCA-------------------
metadata <- ZirFlu$donorSamples %>% full_join(ZirFlu$donorInfo) %>%
  mutate(disease = ifelse(condition == "healthy", "healthy", "cirrhosis")) %>% 
  drop_na() %>%
  mutate(condition = factor(condition, levels = c("healthy", "compensated cirrhosis", "decompensated cirrhosis")))

metabolite_pca <- ZirFlu$metabolite_dat %>% 
  rownames_to_column("probenID") %>% filter(probenID %in% metadata$probenID) %>%
  column_to_rownames("probenID") %>%
  prcomp()
summary(metabolite_pca)$importance[1:3, 1:5]

cowplot::plot_grid(
  get.pca_plot(pca = metabolite_pca, metadat = metadata, groupType = "sex"),
  get.pca_plot(pca = metabolite_pca, metadat = metadata, groupType = "time"),
  nrow = 2
)

get.pca_plot(pca = metabolite_pca, metadat = metadata, groupType = "condition") +
  theme(legend.position = "top")

metadata %>% filter(probenID %in% rownames(ZirFlu$metabolite_dat)) %>% 
  count(season, condition, time)
## metabolite with lm() model between cirrhosis vs. healthy ---------------------------------------------------
library(limma)
library(gplots)

get.limmaRes <- function(metaDat, inputDat) {
  # Aim: identify DE metabolites correct with sex, age, and condition (healthy vs. disease) 
  #by running the linear model (limma package)
  
  # input: metadata - interested participant with sex, age, and conditioin, 
  #inputDat - data with participant in colnames and metabolites in rownames (need to select only overlaped samples)
  
  # output: res - limma output
  ovelapled_samples <- intersect(metaDat$probenID, colnames(inputDat))
  inputDat_temp <- inputDat[, ovelapled_samples]
  
  if (identical(metaDat$probenID, colnames(inputDat_temp)) == TRUE) {
    res <- lmFit(inputDat_temp,
                 design = model.matrix(~ + sex + age + condition, metaDat)) %>% 
      eBayes()
  } else res <- "Error: check input"
  
  return(res)
}

metaboliteDat <- ZirFlu$metabolite_dat %>% 
  rownames_to_column("probenID") %>% 
  filter(probenID %in% ZirFlu$donorSamples$probenID) %>%
  column_to_rownames("probenID") %>% t()

years <- unique(metadata$season)
timepoints <- unique(metadata$time)
res <- list()
for (year in years) {
  for (timepoint in timepoints) {
    res[[paste0(year, "_", timepoint)]] <- get.limmaRes(
      metaDat = metadata %>% filter(season == year, time == timepoint),
      inputDat = metaboliteDat)
  }
}

compareGroupList <- c("conditioncompensated cirrhosis", "conditiondecompensated cirrhosis")
resDE <- res$`2019_T1`$p.value %>% 
  as.data.frame() %>% select(compareGroup) %>% 
  filter(. <0.05) %>%
  rownames_to_column("ionIdx") %>% mutate(ionIdx = as.numeric(ionIdx)) %>%
  left_join(ZirFlu$metabolite_annot %>% select(ionIdx, Formula) %>% distinct())

comparedGroups <- c("conditioncompensated cirrhosis", "conditiondecompensated cirrhosis")
res_padj <- list()
for (i in names(res)) {
  for (comparedGroup in comparedGroups) {
    res_temp <- res[[i]]$p.value %>% as.data.frame() %>% select(comparedGroup) 
    
    outcome_temp <- cbind(res_temp, "p.adj" = p.adjust(unlist(res_temp), method = "fdr")) %>%
      rownames_to_column("ionIdx") %>% mutate(ionIdx = as.numeric(ionIdx)) %>%
      left_join(ZirFlu$metabolite_annot %>% select(ionIdx, Formula) %>% distinct()) %>%
      rename("p.value" := comparedGroup)
    
    res_padj[[i]][[comparedGroup]] <- outcome_temp
  }
}

resDE <- res_padj %>%
  lapply(function(x) x %>% lapply(function(y) y %>% filter(p.value < 0.05)))

resDE_padj <- res_padj %>%
  lapply(function(x) x %>% lapply(function(y) y %>% filter(p.adj < 0.05)))
### Venn diagram --------------------
library(ggvenn)
library(ggVennDiagram)
venn_disease_2019 <- list(
  Baseline_CCvsHealthy = resDE_padj$`2019_T1`$`conditioncompensated cirrhosis`$Formula,
  Baseline_DCvsHealthy = resDE_padj$`2019_T1`$`conditiondecompensated cirrhosis`$Formula,
  Visit3_CCvsHealthy = resDE_padj$`2019_T3`$`conditioncompensated cirrhosis`$Formula,
  Visit3_DCvsHealthy = resDE_padj$`2019_T3`$`conditiondecompensated cirrhosis`$Formula
)

ggvenn(
  venn_disease_2019, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 4) + 
  scale_x_continuous(expand = expansion(mult = .4))

venn_disease_2020 <- list(
  Baseline_CCvsHealthy = resDE_padj$`2020_T1`$`conditioncompensated cirrhosis`$Formula,
  Baseline_DCvsHealthy = resDE_padj$`2020_T1`$`conditiondecompensated cirrhosis`$Formula,
  Visit3_CCvsHealthy = resDE_padj$`2020_T3`$`conditioncompensated cirrhosis`$Formula,
  Visit3_DCvsHealthy = resDE_padj$`2020_T3`$`conditiondecompensated cirrhosis`$Formula
)

ggvenn(
  venn_disease_2020, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 4) + 
  scale_x_continuous(expand = expansion(mult = .4))

# check metabolites

DEmebo_disease_2019 <- c(intersect(venn_disease_2019$Baseline_CCvsHealthy, 
                                  venn_disease_2019$Baseline_DCvsHealthy),
                        intersect(venn_disease_2019$Visit3_CCvsHealthy, 
                                  venn_disease_2019$Visit3_DCvsHealthy))

DEmebo_disease_2020 <- c(intersect(venn_disease_2020$Baseline_CCvsHealthy, 
                                  venn_disease_2020$Baseline_DCvsHealthy),
                        intersect(venn_disease_2020$Visit3_CCvsHealthy, 
                                  venn_disease_2020$Visit3_DCvsHealthy))

## metabolite with lm() model, baseline metabolites to antibody titer at T2, correct to age, gender, and disease ------------------------------------------------
get.lmRes <- function(antibody, metabolites, inputDat) {
  outcome <- data.frame(targetVariable = character(), 
                        independentVariable = character(),
                        estimate = double(), std.error = double(),
                        statistic = double(), p.value = double())
  for (metabolite in metabolites) {
    dat_temp <- inputDat %>% 
      select(all_of(antibody), sex, age, all_of(metabolite), condition) %>%
      rename(antibody = {{antibody}}, metabolite = {{metabolite}})
    
    lm_result <- tidy(lm(antibody ~ sex + age + metabolite + condition, data = dat_temp)) %>%
      select(term, estimate, std.error, statistic, p.value) %>%
      as.data.frame()
    
    outcome[nrow(outcome) + 1, ] <- c(metabolite, lm_result[4, ])
   # outcome[nrow(outcome) + 1, ] <- c(metabolite, lm_result[5, ])
    # outcome[nrow(outcome) + 1, ] <- c(metabolite, lm_result[6, ])
  }
  return(outcome)
}

get.lmRes_padj <- function(res) {
  dat_temp <- p.adjust(res$p.value, method = "fdr")
  lm_padj <- cbind(res, "p.adj" = dat_temp)
  return(lm_padj)
}

# per season
input_metabolite <- list()
input_metabolite[["2019"]] <- ZirFlu$HAItiter_2019 %>% 
  left_join(inputDat_metabolite %>% filter(season == "2019"))
input_metabolite[["2020"]] <- ZirFlu$HAItiter_2020 %>% 
  left_join(inputDat_metabolite %>% filter(season == "2020"))

lmRes <- list()
for (season in c("2019", "2020")) {
  for (antibody in c("H1N1_T2", "H3N2_T2", "Bvictoria_T2", "Byamagata_T2")) {
    lmRes[[paste0(season, "_", antibody)]] <- get.lmRes(
      antibody, 
      metabolites = names(ZirFlu$metabolite_dat), 
      inputDat = input_metabolite[[season]]) %>% 
      get.lmRes_padj() %>%
      mutate(ionIdx = as.numeric(targetVariable)) %>%
      left_join(ZirFlu$metabolite_annot %>% select(ionIdx, Formula) %>% distinct())
  }
}

DE_padj <- lmRes %>%
  lapply(function(x) x %>% filter(p.adj < 0.05)) # only 1 metabolite for 2019_Byamagata_T2

DE_pvalue <- lmRes %>%
  lapply(function(x) x %>% filter(p.value < 0.05))

## Venn diagram --------------------
library(ggvenn)
venn_ab_2019 <- list(
  H1N1 = DE_pvalue$`2019_H1N1_T2`$Formula,
  H3N2 = DE_pvalue$`2019_H3N2_T2`$Formula,
  Bvictoria = DE_pvalue$`2019_Bvictoria_T2`$Formula,
  Byamagata = DE_pvalue$`2019_Byamagata_T2`$Formula
)

ggvenn(
  venn_ab_2019, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 4) + 
  scale_x_continuous(expand = expansion(mult = .4))

venn_ab_2020 <- list(
  H1N1 = DE_pvalue$`2020_H1N1_T2`$Formula,
  H3N2 = DE_pvalue$`2020_H3N2_T2`$Formula,
  Bvictoria = DE_pvalue$`2020_Bvictoria_T2`$Formula,
  Byamagata = DE_pvalue$`2020_Byamagata_T2`$Formula
)

ggvenn(
  venn_ab_2020, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 4) + 
  scale_x_continuous(expand = expansion(mult = .4))

# check metabolites
DEmebo_ab_2019 <- c(intersect(venn_ab_2019$H1N1, 
                              venn_ab_2019$H3N2),
                    intersect(venn_ab_2019$Bvictoria,
                              venn_ab_2019$Byamagata))

DEmebo_ab_2020 <- c(intersect(venn_ab_2020$H1N1, 
                              venn_ab_2020$H3N2),
                      intersect(venn_ab_2020$Bvictoria, 
                                venn_ab_2020$Byamagata))

# heatmap for DE metabolite correlated with antibody titer visit2---------------------
library(ComplexHeatmap)
lmStatistic <- lmRes %>% 
  lapply(function(x) x %>% select(targetVariable, statistic, Formula)) %>%
  imap(~cbind(.x, ab = .y))

lmStatistic <- lmRes %>% 
  imap(~cbind(.x, ab = .y)) %>% purrr::reduce(full_join) %>% 
  mutate(season = substring(ab, 1, 4)) %>%
  mutate(significance = ifelse(p.value < 0.05, TRUE, FALSE))

plotDat_wide <- lmStatistic %>%
  select(statistic, Formula, ab) %>% 
  pivot_wider(names_from = "ab", values_from = statistic) %>%
  column_to_rownames("Formula")

plotDat_wide %>% as.matrix() %>% heatmap(Colv = NA, Rowv = NA)
plotDat_wide %>% as.matrix() %>% Heatmap()
plotDat_wide %>% as.matrix() %>% Heatmap(cluster_columns = FALSE)

lmStatistic %>% 
  ggplot(aes(x = ab, y = Formula)) + 
  geom_tile(aes(fill = statistic)) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Check DE metabolites
DEmebo_2019 <- c(unlist(venn_ab_2019))

# only select metabolites show consistent trend (>75% across strains and years)
consistented_metabolites <- lmStatistic %>% group_by(Formula) %>%
  summarise(postitive = sum(statistic > 0, na.rm = TRUE)) %>%
  full_join(lmStatistic %>% group_by(Formula) %>%
              summarise(negative = sum(statistic < 0, na.rm = TRUE))) %>%
  mutate(select = ifelse(postitive >= 8*0.75 | negative >= 8*0.75, TRUE, NA)) %>%
  filter(select == TRUE)

selected_lmStatistic <- lmStatistic %>% 
  filter(Formula %in% consistented_metabolites$Formula) %>%
  filter(Formula %in% DEmebo_2019)

selected_lmStatistic %>%
  ggplot(aes(x = ab, y = Formula)) + 
  geom_tile(aes(fill = statistic)) +
  geom_text(aes(label = ifelse(significance == TRUE, "*", ""))) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))


# check the boxplot -------------
metaboliteDat <- metadata %>% 
  left_join(ZirFlu$metabolite_dat %>% 
              t() %>% as.data.frame() %>%
              rownames_to_column("ionIdx") %>% mutate(ionIdx = as.numeric(ionIdx)) %>%
              left_join(ZirFlu$metabolite_annot %>% select(ionIdx, Formula) %>% distinct()) %>%
              select(-ionIdx) %>% column_to_rownames("Formula") %>% 
              t() %>% as.data.frame() %>%
              rownames_to_column("probenID"))
metaboliteDat_T1 <- metaboliteDat %>% filter(time == "T1") %>%
  full_join(ZirFlu$HAItiter)

library(ggpubr)
my_comparisons <- list(c("healthy", "compensated cirrhosis"), 
                       c("healthy", "decompensated cirrhosis"), 
                       c("compensated cirrhosis", "decompensated cirrhosis"))
metabolite <- "C5H11NO2"
metabolite <- "C20H27NO3"
metabolite <- "C21H32O6S"
metaboliteDat_T1 %>% 
  ggboxplot(x = "condition", y = metabolite, 
            color = "condition", palette = "jco",
            add = "jitter") +
  stat_compare_means(comparisons = my_comparisons, method = "t.test") +
  facet_wrap("season")

metaboliteDat_T1 %>% 
  ggscatter(x = "H1N1_T2", y = metabolite, palette = "jco",
            add = "reg.line") + stat_cor()
metaboliteDat_T1_T3 <- metaboliteDat %>% filter(time %in% c("T1", "T3")) %>%
  full_join(ZirFlu$HAItiter)

metabolite <- "C21H32O6S"
metaboliteDat_T1_T3 %>%
  ggboxplot(x = "condition", y = metabolite, 
            color = "condition", palette = "jco",
            add = "jitter") +
  stat_compare_means(comparisons = my_comparisons, method = "t.test") +
  facet_grid(season~time)

plot_temp <- cowplot::plot_grid(
  metaboliteDat_T1_T3 %>% filter(season == "2019") %>%
    ggboxplot(x = "condition", y = "C21H32O6S", 
              color = "condition", palette = "jco",
              add = "jitter") +
    stat_compare_means(comparisons = my_comparisons, method = "t.test") +
    facet_wrap("time"),
  metaboliteDat_T1_T3 %>% filter(season == "2019") %>%
    ggboxplot(x = "condition", y = "C20H27NO3", 
              color = "condition", palette = "jco",
              add = "jitter") +
    stat_compare_means(comparisons = my_comparisons, method = "t.test") +
    facet_wrap("time"),
  metaboliteDat_T1_T3 %>% filter(season == "2019") %>%
    ggboxplot(x = "condition", y ="C39H77O10P", 
              color = "condition", palette = "jco",
              add = "jitter") +
    stat_compare_means(comparisons = my_comparisons, method = "t.test") +
    facet_wrap("time"),
  nrow = 1
)
pdf(file = "20230428_pickedMebo_T1vsT3.pdf", width = 16, height = 4, onefile = TRUE)
plot(plot_temp)
dev.off()

metaboliteDat_T1_T3 %>% filter(season == "2019") %>%
  ggboxplot(x = "condition", y = metabolite, 
            color = "condition", palette = "jco",
            add = "jitter") +
  stat_compare_means(comparisons = my_comparisons, method = "t.test") +
  facet_wrap("time")


pdf(file = "20230427_metaboliteT1_abT2.pdf", width = 16, height = 8, onefile = TRUE)
for(metabolite in unique(selected_lmStatistic$Formula)) {
  plot_temp <- cowplot::plot_grid(
    # season 2019
    cowplot::plot_grid(
      metaboliteDat_T1 %>% filter(season == "2019") %>%
        ggscatter(x = "H1N1_T2", y = metabolite, color = "condition", palette = "jco",
                  add = "reg.line") + stat_cor(aes(color = condition)),
      metaboliteDat_T1 %>% filter(season == "2019") %>%
        ggscatter(x = "H3N2_T2", y = metabolite, color = "condition", palette = "jco",
                  add = "reg.line") + stat_cor(aes(color = condition)),
      metaboliteDat_T1 %>% filter(season == "2019") %>%
        ggscatter(x = "Bvictoria_T2", y = metabolite, color = "condition", palette = "jco",
                  add = "reg.line") + stat_cor(aes(color = condition)),
      metaboliteDat_T1 %>% filter(season == "2019") %>%
        ggscatter(x = "Byamagata_T2", y = metabolite, color = "condition", palette = "jco",
                  add = "reg.line") + stat_cor(aes(color = condition)),
      nrow = 1
    ),
    # season 2020
    cowplot::plot_grid(
      metaboliteDat_T1 %>% filter(season == "2020") %>%
        ggscatter(x = "H1N1_T2", y = metabolite, color = "condition", palette = "jco",
                  add = "reg.line") + stat_cor(aes(color = condition)),
      metaboliteDat_T1 %>% filter(season == "2020") %>%
        ggscatter(x = "H3N2_T2", y = metabolite, color = "condition", palette = "jco",
                  add = "reg.line") + stat_cor(aes(color = condition)),
      metaboliteDat_T1 %>% filter(season == "2020") %>%
        ggscatter(x = "Bvictoria_T2", y = metabolite, color = "condition", palette = "jco",
                  add = "reg.line") + stat_cor(aes(color = condition)),
      metaboliteDat_T1 %>% filter(season == "2020") %>%
        ggscatter(x = "Byamagata_T2", y = metabolite, color = "condition", palette = "jco",
                  add = "reg.line") + stat_cor(aes(color = condition)),
      nrow = 1
    ),
    nrow = 2
  )
  plot(plot_temp )
}
dev.off()

# pathway analysis at cpdb ----------------------------------
library("KEGGREST")

kegg_annot <- c() # KEGG database: using formula to annotate metabolite 
for (formula in unique(selected_lmStatistic$Formula)) {
  annot_temp <- keggFind("compound", formula, "formula")
  exact_formula <- annot_temp[which(annot_temp == formula)]
  kegg_annot <- c(kegg_annot, exact_formula)
} # 71 formulas (183 Ids)

kegg_outcomes <- kegg_annot %>% enframe() %>% distinct() %>%
  rename("cpdId" = "name", "Formula" = "value") %>% 
  separate(col = "cpdId", into = c("cpd", "keggIds"), sep = ":") # 31 formulas were able to annotate

write.table(kegg_outcomes$keggIds, 
     file = "temp/20230428_keggID_DE_meboT1_abT2.txt",
     col.names = FALSE, row.names = FALSE, quote = FALSE)
