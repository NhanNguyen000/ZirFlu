rm(list = ls())

library(tidyverse)
library(limma)
library(gplots)
library(ggvenn)
library(ggVennDiagram)
library(broom)
# functions used in the analysis -----------------------------------------------

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


# load the data list object------------------------------------------------------
load("data/ZirFlu.RData")

# lm() model ------------------------------------------------------------------
# Aim: calculate the association  between baseline metabolite to antibody titer at Visit 2, 
# after correcting to age, gender, and disease

metadata <- ZirFlu$donorSamples %>% 
  full_join(ZirFlu$donorInfo) %>%
  drop_na() %>%
  mutate(condition = factor(condition, 
                            levels = c("Healthy", "Compensated cirrhosis", "Decompensated cirrhosis")))


meboDat_baseline <- ZirFlu$metabolite_dat %>% 
  rownames_to_column("probenID") %>% 
  right_join(ZirFlu$donorSamples %>% filter(time == "Baseline"))

input_metabolite <- ZirFlu$HAItiter %>% 
  full_join(ZirFlu$donorInfo) %>% 
  inner_join(meboDat_baseline)

years <- unique(metadata$season)
ab_levels <- c("H1N1_Visit2", "H3N2_Visit2", "Bvictoria_Visit2", "Byamagata_Visit2")

lmRes <- list()
for (year in years) {
  for (ab_level in ab_levels){
    lmRes[[paste0(year, "_", ab_level)]] <- get.lmRes(
      antibody = ab_level, 
      metabolites = names(ZirFlu$metabolite_dat), 
      inputDat = input_metabolite %>% filter(season == year)) %>% 
      get.lmRes_padj() %>%
      mutate(ionIdx = as.numeric(targetVariable)) %>%
      left_join(ZirFlu$metabolite_annot %>% select(ionIdx, Formula) %>% distinct())
  }
}

## save lm() model outcome -------------------------------------------------------------------
save(lmRes, file = "data/lmRes_meboAbVisit2.RData")

## check sig. metabolites  -------------------------------------------------------------------
DE_padj <- lmRes %>%
  lapply(function(x) x %>% filter(p.adj < 0.05)) # only 1 metabolite for 2019_Byamagata_Visit2

DE_pvalue <- lmRes %>%
  lapply(function(x) x %>% filter(p.value < 0.05))

# venn diagram for overlapped metabolites --------------------------------------

## season 2019 ------------------------------------------------------------------
venn_meboAbVisit2_2019 <- list(
  Bvictoria = DE_pvalue$`2019_Bvictoria_Visit2`$Formula,
  H1N1 = DE_pvalue$`2019_H1N1_Visit2`$Formula,
  H3N2 = DE_pvalue$`2019_H3N2_Visit2`$Formula,
  Byamagata = DE_pvalue$`2019_Byamagata_Visit2`$Formula
)

vennPlot_meboAbVisit2_2019 <- ggvenn(
  venn_meboAbVisit2_2019, 
  fill_color = c("#FBE2D1", "#D3C3BE", "#FFF9E5", "#E7E7E7"),
  stroke_size = 0.3, set_name_size = 5, text_size = 6) + 
  scale_x_continuous(expand = expansion(mult = .2))

# Save the plot as SVG
vennPlot_meboAbVisit2_2019
ggsave("output/vennPlot_meboAbVisit2_seasons2019.svg", 
       vennPlot_meboAbVisit2_2019, device = "svg")
ggsave("output/vennPlot_meboAbVisit2_seasons2019.png", 
       vennPlot_meboAbVisit2_2019, device = "png")


## season 2020 ------------------------------------------------------------------
venn_meboAbVisit2_2020 <- list(
  Bvictoria = DE_pvalue$`2020_Bvictoria_Visit2`$Formula,
  H1N1 = DE_pvalue$`2020_H1N1_Visit2`$Formula,
  H3N2 = DE_pvalue$`2020_H3N2_Visit2`$Formula,
  Byamagata = DE_pvalue$`2020_Byamagata_Visit2`$Formula
)

vennPlot_meboAbVisit2_2020 <- ggvenn(
  venn_meboAbVisit2_2020, 
  fill_color = c("#FBE2D1", "#D3C3BE", "#FFF9E5", "#E7E7E7"),
  stroke_size = 0.3, set_name_size = 5, text_size = 6) + 
  scale_x_continuous(expand = expansion(mult = .2))

# Save the plot as SVG
vennPlot_meboAbVisit2_2020
ggsave("output/vennPlot_meboAbVisit2_seasons2020.svg", 
       vennPlot_meboAbVisit2_2020, device = "svg")
ggsave("output/vennPlot_meboAbVisit2_seasons2020.png", 
       vennPlot_meboAbVisit2_2020, device = "png")
## save DE metabolites ----------------------------------------------------------
save(venn_meboAbVisit2_2019, venn_meboAbVisit2_2020,
     file = "data/DEmebo_abVisit2.RData")

# check metabolites -------------------------------------------------------------
# option 1:
DEmebo_abVisit2_2019 <- intersect(intersect(venn_meboAbVisit2_2019$H1N1, 
                                    venn_meboAbVisit2_2019$H3N2),
                          intersect(venn_meboAbVisit2_2019$Bvictoria, 
                                    venn_meboAbVisit2_2019$Byamagata))

DEmebo_abVisit2_2020 <- intersect(intersect(venn_meboAbVisit2_2020$H1N1, 
                                    venn_meboAbVisit2_2020$H3N2),
                          intersect(venn_meboAbVisit2_2020$Bvictoria, 
                                    venn_meboAbVisit2_2020$Byamagata))

# option 2: 
DEmebo_abVisit2_2019 <- c(intersect(venn_meboAbVisit2_2019$H1N1, 
                                    venn_meboAbVisit2_2019$H3N2),
                         intersect(venn_meboAbVisit2_2019$Bvictoria, 
                                   venn_meboAbVisit2_2019$Byamagata))

DEmebo_abVisit2_2020 <- c(intersect(venn_meboAbVisit2_2020$H1N1, 
                                    venn_meboAbVisit2_2020$H3N2),
                          intersect(venn_meboAbVisit2_2020$Bvictoria, 
                                    venn_meboAbVisit2_2020$Byamagata))

