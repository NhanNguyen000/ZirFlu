library(tidyverse)
library(lme4)
library(lmerTest)
library(broom.mixed)
library(ggpubr)

# protein -----------------
inputDat <- ZirFlu$donorInfo %>% full_join(ZirFlu$donorSamples) %>%
  select(patientID, probenID, sex, age, condition, category, time) %>%
  full_join(ZirFlu$HAItiter_Feb2022 %>% select(-condition)) %>%
  right_join(ZirFlu$protein_dat %>% rownames_to_column("probenID") %>% mutate(probenID = as.numeric(probenID)))


lm_results <- list()
var_associDisease <- list()
var_associAbTiter <- list()
var_onlyAbTiter <- list()
timeList <- c("T1", "T3", "T4")
abTiters <- c("H1N1_titerFC", "H3N2_titerFC", 
              "Bvictoria.Maryland_titerFC", "Byamagata.Phuket_titerFC")
varList <- names(ZirFlu$protein_dat)
for (timepoint in timeList) {
  dat_time <- inputDat %>% filter(time %in% timepoint)
  
  for (name in abTiters) {
    resName <- paste0(timepoint, "_", name)
    
    lm_results[[resName]] <- get.lmTest_correction(variableList = varList,
                                                abTiter = name, inputDat = dat_time)
    var_associDisease[[resName]] <- lm_results[[resName]] %>% filter(p.value < 0.05) %>%
      filter(independentVariable == "diseasehealthy")
    
    var_associAbTiter[[resName]] <- lm_results[[resName]] %>% filter(p.value < 0.05) %>%
      filter(independentVariable == "abTiter")
    
    var_onlyAbTiter[[resName]] <- get.var_onlyAbTiter(var_associAbTiter, 
                                                      var_associDisease,
                                                      abTiter = resName)
  }
}

save(lm_results, var_associAbTiter, var_associDisease, 
     var_onlyAbTiter, file = "temp/lmResult.RData")

# test padjust value:
lm_padj <- list()
var_associDisease_padj <- list()
var_associAbTiter_padj <- list()
for (i in names(lm_results)) {
  dat_temp <- p.adjust(lm_results[[i]]$p.value, method = "fdr")
  lm_padj[[i]] <- cbind(lm_results[[i]], "p.adj" = dat_temp)
  
  var_associDisease_padj[[i]] <- lm_padj[[i]] %>% filter(p.value < 0.05) %>%
    filter(independentVariable == "diseasehealthy")
  
  var_associAbTiter_padj[[i]] <- lm_padj[[i]] %>% filter(p.value < 0.05) %>%
    filter(independentVariable == "abTiter")
}

## venn diagram -------------------
load("temp/lmResult.RData")
library(ggVennDiagram)

get.vennDat <- function(inputDat, time) {
  dat_temp <- inputDat %>% lapply(function(x) x%>% select(targetVariable))
  outcome <- dat_temp[grep(time, names(dat_temp))] %>% unlist(recursive = FALSE)
  names(outcome) <- substring(names(outcome), 1, 7)
  return(outcome)
}

cowplot::plot_grid(
  ggVennDiagram(get.vennDat(var_associDisease, "T1"),
                label_alpha = 0) + scale_fill_gradient(low="white",high = "blue"),
  ggVennDiagram(get.vennDat(var_associAbTiter, "T1"),
                label_alpha = 0) + scale_fill_gradient(low="white",high = "blue"),
  labels = c("DEPs disease", "DEPs abTiter"), nrow = 2)

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

# heatmap --------------------
get.lmStatistic <- function(inputDat) {
  outcome <- inputDat %>% 
    lapply(function(x) x %>% select(targetVariable, statistic)) %>% 
    bind_rows(.id = "groups") %>% mutate(groups = substring(groups, 1, 7)) %>%
    pivot_wider(names_from = groups, values_from = statistic) %>% 
    column_to_rownames("targetVariable")
  return(outcome)
}

varDisease <- get.lmStatistic(var_associDisease) %>% mutate_if(is.numeric, as.logical)
varAbTiter <- get.lmStatistic(var_associAbTiter) %>% mutate_if(is.numeric, as.logical)
varlmAbTiter <- get.lmStatistic(lm_results %>% 
                                      lapply(function(x) x%>% filter(independentVariable == "abTiter"))) %>%
  as.matrix()

# heatmap
heatmap(varlmAbTiter, Colv = NA, Rowv = NA)
varlmAbTiter[which(rownames(varlmAbTiter) %in% rownames(varAbTiter)),] %>%
  heatmap(Colv = NA, Rowv = NA)



heatmapDat <- varlmAbTiter[which(rownames(varlmAbTiter) %in% unique(unlist(var_onlyAbTiter))),] %>% 
  as.data.frame() %>% 
  rownames_to_column("OlinkID") %>% 
  pivot_longer(cols = starts_with("T"), 
               names_to = "time", values_to = "statistic")

heatmapDat %>% ggplot(aes(x = time, y = OlinkID)) + 
  geom_tile(aes(fill = statistic)) +
  scale_fill_continuous(low = "violetred", high = "aquamarine") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

outcome <- as.data.frame(matrix(nrow=0, ncol = 4))
for (name in names(var_onlyAbTiter)) {
  
  dat_temp <- heatmapDat %>% filter(time == substring(name, 1, 7)) %>% 
    mutate(significance = ifelse(OlinkID %in% var_onlyAbTiter[[name]], TRUE, NA))
  
  outcome <- rbind(outcome, dat_temp)
}
all.equal(heatmapDat %>% arrange(time), 
           outcome %>% arrange(time) %>% select(1:3)) # TRUE = same data
outcome %>% left_join(ZirFlu$var_annot) %>% 
  ggplot(aes(x = time, y = Assay)) + 
  geom_tile(aes(fill = statistic)) +
  geom_text(aes(label = ifelse(significance == TRUE, "*", ""))) +
  scale_fill_continuous(low = "violetred", high = "aquamarine") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# check complex heatmap package

# metabolite -----------------
inputDat <- ZirFlu$donorInfo %>% full_join(ZirFlu$donorSamples) %>%
  select(patientID, probenID, sex, age, condition, category, time) %>%
  full_join(ZirFlu$HAItiter_Feb2022 %>% select(-condition)) %>%
  right_join(ZirFlu$metabolite_dat %>% rownames_to_column("probenID") %>% mutate(probenID = as.numeric(probenID)))


lm_results <- list()
var_associDisease <- list()
var_associAbTiter <- list()
var_onlyAbTiter <- list()
timeList <- c("T1", "T3", "T4")
abTiters <- c("H1N1_titerFC", "H3N2_titerFC", 
              "Bvictoria.Maryland_titerFC", "Byamagata.Phuket_titerFC")
varList <- names(ZirFlu$metabolite_dat)
for (timepoint in timeList) {
  dat_time <- inputDat %>% filter(time %in% timepoint)
  
  for (name in abTiters) {
    resName <- paste0(timepoint, "_", name)
    
    lm_results[[resName]] <- get.lmTest_correction(variableList = varList,
                                                   abTiter = name, inputDat = dat_time)
    var_associDisease[[resName]] <- lm_results[[resName]] %>% filter(p.value < 0.05) %>%
      filter(independentVariable == "diseasehealthy")
    
    var_associAbTiter[[resName]] <- lm_results[[resName]] %>% filter(p.value < 0.05) %>%
      filter(independentVariable == "abTiter")
    
    var_onlyAbTiter[[resName]] <- get.var_onlyAbTiter(var_associAbTiter, 
                                                      var_associDisease,
                                                      abTiter = resName)
  }
}

save(lm_results, var_associAbTiter, var_associDisease, 
     var_onlyAbTiter, file = "temp/lmResult_metabolitesNoFilter.RData")

# test padjust value:
lm_padj <- list()
var_associDisease_padj <- list()
var_associAbTiter_padj <- list()
for (i in names(lm_results)) {
  dat_temp <- p.adjust(lm_results[[i]]$p.value, method = "fdr")
  lm_padj[[i]] <- cbind(lm_results[[i]], "p.adj" = dat_temp)
  
  var_associDisease_padj[[i]] <- lm_padj[[i]] %>% filter(p.value < 0.05) %>%
    filter(independentVariable == "diseasehealthy")
  
  var_associAbTiter_padj[[i]] <- lm_padj[[i]] %>% filter(p.value < 0.05) %>%
    filter(independentVariable == "abTiter")
}

## venn diagram -------------------
load("temp/lmResult_metabolitesNoFilter.RData")
library(ggVennDiagram)

get.vennDat <- function(inputDat, time) {
  dat_temp <- inputDat %>% lapply(function(x) x%>% select(targetVariable))
  outcome <- dat_temp[grep(time, names(dat_temp))] %>% unlist(recursive = FALSE)
  names(outcome) <- substring(names(outcome), 1, 7)
  return(outcome)
}

cowplot::plot_grid(
  ggVennDiagram(get.vennDat(var_associDisease, "T1"),
                label_alpha = 0) + scale_fill_gradient(low="white",high = "blue"),
  ggVennDiagram(get.vennDat(var_associAbTiter, "T1"),
                label_alpha = 0) + scale_fill_gradient(low="white",high = "blue"),
  labels = c("DEPs disease", "DEPs abTiter"), nrow = 2)

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

# mlic and citric acid expression
citric - C6H8O7
malate - C4H6O5 - ionId 132

## T1 protein plot -------------------
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
