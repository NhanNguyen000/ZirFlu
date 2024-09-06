rm(list = ls())

library(tidyverse)
library(limma)
library(gplots)
library(ggvenn)
library(ggVennDiagram)
library(openxlsx)
# functions used in the analysis -----------------------------------------------

get.limmaRes <- function(metaDat, inputDat) {
  # Aim: estimate the relationship between metabolites level and sex, age, health condition (healthy vs. disease) 
  #by running the linear model (limma package)
  
  # input: metadata - interested participant with sex, age, and condition, 
  # inputDat - data with participant in colnames and metabolites in rownames (need to select only overlaped samples)
  
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

# load the data list object------------------------------------------------------
load("processedData/ZirFlu.RData")

# lm() model ------------------------------------------------------------------
# Aim: calculate the association between metabolite and health conditions (cirrhosis vs. healthy), 
# after correcting to age and gender

metadata <- ZirFlu$donorSamples %>% 
  full_join(ZirFlu$donorInfo) %>%
  drop_na() %>%
  mutate(condition = factor(condition, 
        levels = c("Healthy", "Compensated cirrhosis", "Decompensated cirrhosis")))


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


## save lm() model outcome -------------------------------------------------------------------
save(res, file = "processedData/lmRes_disease.RData")

## check sig. metabolites  -------------------------------------------------------------------
comparedGroups <- c("conditionCompensated cirrhosis", "conditionDecompensated cirrhosis")
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

## save as excel files
write.xlsx(list(Baseline_2019_compCirrhosis = res_padj$`2019_Baseline`$`conditionCompensated cirrhosis`,
                Baseline_2019_decompCirrhosis = res_padj$`2019_Baseline`$`conditionDecompensated cirrhosis`,
                Visit1_2019_compCirrhosis = res_padj$`2019_Visit1`$`conditionCompensated cirrhosis`,
                Visit1_2019_decompCirrhosis = res_padj$`2019_Visit1`$`conditionDecompensated cirrhosis`,
                Visit2_2019_compCirrhosis = res_padj$`2019_Visit2`$`conditionCompensated cirrhosis`,
                Visit2_2019_decompCirrhosis = res_padj$`2019_Visit2`$`conditionDecompensated cirrhosis`,
                Baseline_2020_compCirrhosis = res_padj$`2020_Baseline`$`conditionCompensated cirrhosis`,
                Baseline_2020_decompCirrhosis = res_padj$`2020_Baseline`$`conditionDecompensated cirrhosis`,
                Visit1_2020_compCirrhosis = res_padj$`2020_Visit1`$`conditionCompensated cirrhosis`,
                Visit1_2020_decompCirrhosis = res_padj$`2020_Visit1`$`conditionDecompensated cirrhosis`,
                Visit2_2020_compCirrhosis =  res_padj$`2020_Visit2`$`conditionCompensated cirrhosis`,
                Visit2_2020_decompCirrhosis = res_padj$`2020_Visit2`$`conditionDecompensated cirrhosis`), 
           "output/suppTables_metabolite_associatedDiseaseCondition.xlsx", rowNames = FALSE)

# Venn diagram for overlapped metabolites --------------------------------------

## season 2019 ------------------------------------------------------------------
venn_disease_2019 <- list(
  Baseline_CCvsHealthy = resDE_padj$`2019_Baseline`$`conditionCompensated cirrhosis`$Formula,
  Baseline_DCvsHealthy = resDE_padj$`2019_Baseline`$`conditionDecompensated cirrhosis`$Formula,
  Visit2_CCvsHealthy = resDE_padj$`2019_Visit2`$`conditionCompensated cirrhosis`$Formula,
  Visit2_DCvsHealthy = resDE_padj$`2019_Visit2`$`conditionDecompensated cirrhosis`$Formula
)

vennPlot_disease_2019 <- ggvenn(
  venn_disease_2019, 
  fill_color = c("#FBE2D1", "#D3C3BE", "#FFF9E5", "#E7E7E7"),
  stroke_size = 0.3, set_name_size = 5, text_size = 6) + 
  scale_x_continuous(expand = expansion(mult = .2))


# Save the plot as SVG
vennPlot_disease_2019
ggsave("output/vennPlot_disease_season2019.svg", 
       vennPlot_disease_2019, device = "svg")

ggsave("output/vennPlot_disease_season2019.png", 
       vennPlot_disease_2019, device = "png")

## season 2020 ------------------------------------------------------------------
venn_disease_2020 <- list(
  Baseline_CCvsHealthy = resDE_padj$`2020_Baseline`$`conditionCompensated cirrhosis`$Formula,
  Baseline_DCvsHealthy = resDE_padj$`2020_Baseline`$`conditionDecompensated cirrhosis`$Formula,
  Visit2_CCvsHealthy = resDE_padj$`2020_Visit2`$`conditionCompensated cirrhosis`$Formula,
  Visit2_DCvsHealthy = resDE_padj$`2020_Visit2`$`conditionDecompensated cirrhosis`$Formula
)

vennPlot_disease_2020 <- ggvenn(
  venn_disease_2020, 
  fill_color = c("#FBE2D1", "#D3C3BE", "#FFF9E5", "#E7E7E7"),
  stroke_size = 0.3, set_name_size = 5, text_size = 6) + 
  scale_x_continuous(expand = expansion(mult = .2))

# Save the plot as SVG
vennPlot_disease_2020
ggsave("output/vennPlot_disease_season2020.svg", 
       vennPlot_disease_2020, device = "svg")
ggsave("output/vennPlot_disease_season2020.png", 
       vennPlot_disease_2020, device = "png")

## save DE metabolites ----------------------------------------------------------
save(venn_disease_2019, venn_disease_2020,
     file = "processedData/DEmebo_disease.RData")

# check metabolites ---------------------------------------------------------------
# option 1: 
DEmebo_disease_2019 <- intersect(intersect(venn_disease_2019$Baseline_CCvsHealthy, 
                                   venn_disease_2019$Baseline_DCvsHealthy),
                         intersect(venn_disease_2019$Visit2_CCvsHealthy, 
                                   venn_disease_2019$Visit2_DCvsHealthy))

DEmebo_disease_2020 <- intersect(intersect(venn_disease_2020$Baseline_CCvsHealthy, 
                                   venn_disease_2020$Baseline_DCvsHealthy),
                         intersect(venn_disease_2020$Visit2_CCvsHealthy, 
                                   venn_disease_2020$Visit2_DCvsHealthy))

intersect(DEmebo_disease_2019, DEmebo_disease_2020) # "C19H28O5S"  "C19H30O5S"  "C26H45NO6S" "C26H45NO7S"

# option 2: 
DEmebo_disease_2019 <- c(intersect(venn_disease_2019$Baseline_CCvsHealthy, 
                                   venn_disease_2019$Baseline_DCvsHealthy),
                         intersect(venn_disease_2019$Visit3_CCvsHealthy, 
                                   venn_disease_2019$Visit3_DCvsHealthy))

DEmebo_disease_2020 <- c(intersect(venn_disease_2020$Baseline_CCvsHealthy, 
                                   venn_disease_2020$Baseline_DCvsHealthy),
                         intersect(venn_disease_2020$Visit3_CCvsHealthy, 
                                   venn_disease_2020$Visit3_DCvsHealthy))
