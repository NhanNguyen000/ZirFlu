library(tidyverse)
library(rstatix)
library(ggpubr)
library(stringr)
library(caret)
# overview ---------------------------------------------------------------------
ZirFlu$HAItiter %>% count(season)
ZirFlu$HAItiter %>% filter(season == "2019") %>% count(condition)
ZirFlu$HAItiter %>% filter(season == "2020") %>% count(condition)

# NOTE:
# The participant Z-62 season 2019 were excluded because no serum samples were available.
# The participant Z-01-99-069 season 2020 has azathioprine in his medication 
# and were excluded from the analysis.

HAItiter <- list()
for (i in c("abFC", "T1", "T2", "T3")) {
  HAItiter[[i]] <- ZirFlu$HAItiter %>% 
    select(patientID, season, condition,  disease, matches(i)) %>%
    pivot_longer(!c("patientID", "season", "condition", "disease"), 
                 names_to = "strains", values_to = "log2_value")
}

# Compare 2 cirrhosis groups: decompensated vs. compensated cirrhosis ----------
## boxplot - HAI titer or ab titer ---------------------------------------------

### HAI titer
datPlot <- HAItiter$abFC %>% filter(disease == "cirrhosis") # HAI_abFC (T2 vs. T1)
datPlot <- HAItiter$T1 %>% filter(disease == "cirrhosis") # HAI titer at T1
datPlot <- HAItiter$T2 %>% filter(disease == "cirrhosis") # HAI titer at T2
datPlot <- HAItiter$T3 %>% filter(disease == "cirrhosis") # HAI titer at T3

### ab Titer
datPlot <- HAItiter$T2 %>% filter(disease == "cirrhosis") %>% 
  filter(str_detect(strains, "Byamagata"))
  
bxp <- datPlot %>%
  ggboxplot(x = "season", y = "log2_value", fill = "condition",
            bxp.errorbar = TRUE, palette = "npg",
            outlier.shape = 1, outlier.size = 4, outlier.color = "grey")

stat.test_condition <- datPlot %>% group_by(season) %>%
  t_test(log2_value ~ condition)%>% add_xy_position(x = "season", dodge = 0.8)

stat.test_condition_strain <- datPlot %>% group_by(season, strains) %>%
  t_test(log2_value ~ condition)%>% add_xy_position(x = "season", dodge = 0.8)

stat.test_2seasons <- datPlot %>%
  t_test(log2_value ~ condition)%>% add_xy_position(x = "season", dodge = 0.8)

stat.test_2seasons_strain <- datPlot %>% group_by(strains) %>%
  t_test(log2_value ~ condition)%>% add_xy_position(x = "season", dodge = 0.8)

bxp + 
  stat_pvalue_manual(
    stat.test_condition, label = "p", tip.length = 0.01, bracket.nudge.y = 1)  + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) + ylim(0, 15) + 
  ylab("log2(ab_FC)")
#ylab("log2(HAI_T1)")
#ylab("log2(HAI_T2)")
#ylab("log2(HAI_T3)")
#ylab("log2(HAI_T2_Bvictoria)")

## meta-analysis (?, check previous version) -----------------------------------

# make PCA for protein and metabolite data -------------------------------------------
## HAI ------------------------------------------------------------
metadata <- ZirFlu$HAItiter
# metadata <- ZirFlu$HAItiter %>% filter(season == "2019") # per season
HAItiter_raw <- ZirFlu$HAItiter %>% 
  select(-c(patientID, season, vaccine_response, category, condition, disease))

HAItiter_impute <- predict(preProcess(x=HAItiter_raw, method = "knnImpute"), 
                           HAItiter_raw)
HAItiter_impute %>% 
  prcomp() %>% get.pca_plot(metadata, "condition")

## protein ------------------------------------------------------------
metadata <- ZirFlu$donorSamples %>% full_join(ZirFlu$donorInfo) %>% 
  filter(patientID %in% ZirFlu$HAItiter$patientID) %>%
  filter(probenID %in% rownames(ZirFlu$protein_dat)) %>% 
  na.omit() %>% as.data.frame() # %>% filter(season == "2019") # per season

protein_raw <- ZirFlu$protein_dat %>% 
  rownames_to_column(var = "probenID") %>% 
  filter(probenID %in% metadata$probenID) %>%
  column_to_rownames("probenID")

protein_impute <- 
  predict(preProcess(x=protein_raw, method = "knnImpute"), protein_raw)

# all protein at all time points
protein_impute %>% prcomp() %>% get.pca_plot(metadata, "condition")
protein_impute %>% prcomp() %>% get.pca_plot(metadata, "disease")

# protein at T1 - baseline
metadata_T1 <- metadata %>% filter(time == "T1")
protein_impute_T1 <- protein_impute %>% rownames_to_column(var = "probenID") %>% 
  filter(probenID %in% metadata_T1$probenID) %>% column_to_rownames("probenID")
protein_impute_T1 %>% prcomp() %>% get.pca_plot(metadata_T1, "disease")

## metabolite ------------------------------------------------------------
metabolite_raw <- ZirFlu$metabolite_dat %>% rownames_to_column(var = "probenID") %>%
  filter(probenID %in% metadata$probenID) %>% column_to_rownames("probenID")

# all metabolites at all time points
metabolite_raw %>% prcomp() %>% get.pca_plot(metadata, "condition")
metabolite_raw %>% prcomp() %>% get.pca_plot(metadata, "disease")

# metabolites at T1 - baseline
metabolite_raw_T1 <- metabolite_raw %>% rownames_to_column(var = "probenID") %>% 
  filter(probenID %in% metadata_T1$probenID) %>% column_to_rownames("probenID")
metabolite_raw_T1 %>% 
  prcomp() %>% get.pca_plot(metadata_T1, "disease")
