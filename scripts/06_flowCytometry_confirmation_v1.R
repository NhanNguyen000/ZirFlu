library(readxl)
library(tidyverse)
library(openxlsx)

# functions ------------------------------
get.celltype_metadata <- function(rawDat) {
  metadata <- rawDat[1, -c(1:2)] %>% 
    t() %>% as.data.frame() %>%
    rownames_to_column(var="cell_type") %>%
    rename("subcell_type" = "V1")
  
  metadata$cell_type <- c(rep("CD4+ T cells", 16), 
                          rep("CD8+ T cells", 9),
                          rep("Cytokine secreting T cells", 10), 
                          rep("B cells", 16))
  return(metadata)
}

get.celltype_data <- function(rawDat) {
  celltype_data <- rawDat
  colnames(celltype_data) <- c( "condition", celltype_data[1, -1])
  
  celltype_data <- celltype_data %>% 
    slice(-1) %>% drop_na("Donor ID") %>%
    fill(condition) %>% as.data.frame()
  
  celltype_data[, 3:53] <- lapply(celltype_data[, 3:53], as.numeric)
  
  return(celltype_data)
}
#  available samples ------------------------------
flowCyto_metadata <- list()
flowCyto_data <- list()

sheet_names <- c("Baseline ex vivo", "Baseline re-stimulated", 
                 "Visit 2 ex vivo", "Visit 2 re-stimulated")

for (sheet_name in sheet_names) {
  dat_temp <- read_xlsx("../flowCytometry/20230105_ZirFlu2019-2020_FlowCytometry_Valerie&Nhan.xlsx",
                        sheet = sheet_name)
  
  flowCyto_metadata[[sheet_name]] <- get.celltype_metadata(dat_temp)
  flowCyto_data[[sheet_name]] <- get.celltype_data(dat_temp) %>% 
    mutate(condition2 = factor(ifelse(condition == "Healthy", "healthy", "cirrhosis")))
}

library(ggplot2)
library(ggpubr)

get.ggboxplot <- function(dat, x_col, y_col) {
  ggboxplot(dat, x = x_col, y = y_col,
            color = "condition2", 
            paletter = "jco", add = "jitter") +
    stat_compare_means(method = "t.test")
}

plot_list <- list()
for(celltype in flowCyto_metadata$`Baseline ex vivo`$subcell_type) {
  plot_list[[celltype]] <- cowplot::plot_grid(
    get.ggboxplot(flowCyto_data$`Baseline ex vivo`,
                  x_col = "condition2", y_col = celltype),
    get.ggboxplot(flowCyto_data$`Baseline re-stimulated`,
                  x_col = "condition2", y_col = celltype),
    get.ggboxplot(flowCyto_data$`Visit 2 ex vivo`,
                  x_col = "condition2", y_col = celltype),
    get.ggboxplot(flowCyto_data$`Visit 2 re-stimulated`,
                  x_col = "condition2", y_col = celltype),
    nrow = 1
  )
}

pdf(file = "temp/20230206_cytoflowmetry_boxplot_2019.pdf", 
    width = 12, height = 4 , onefile = TRUE)
for (celltype in names(plot_list)) {
  plot(plot_list[[celltype]])
}
dev.off()

plot_list$Th17
plot_list$Th1
plot_list$Th2
plot_list$`CD4+ Central Memory`
plot_list$`CD4+ Effector`
plot_list$`CD4+ IL-10 Tregs`
plot_list$`CD8+ Effector Memory`
plot_list$`CD8+ Effector`
plot_list$`CD8+ CD28+`
plot_list$`CD8+ CD38+`
plot_list$`CD4+ IL-10+`
plot_list$`CD8+ IFNy+`
plot_list$`CD8+ IL-10+`
plot_list$`B cells`
plot_list$`CD80 activated`
plot_list$`CD86 activated`
plot_list$Naive
plot_list$Memory
plot_list$`IgA switched`
plot_list$`IgG switched`
plot_list$`IgM non-switched`
plot_list$Bregs
plot_list$`IL35 Bregs`
plot_list$Transitional
plot_list$`Plasma cells`
plot_list$Plasmablasts
# check if protein & cell type associated -----------------
# get sample at baseline
sample_T1 <- ZirFlu$donorSamples %>% filter(time == "T1" & season == "2019")

# proteins
selected_proteins <- ZirFlu$protein_annot %>% 
  filter(Assay %in% c("TNFSF10", "CXCL8", "IL6", "IL17C", "IL17D", "TNFSF12"))

protein_dat <- ZirFlu$protein_dat %>% select(selected_proteins$OlinkID) %>% 
  rownames_to_column("probenID") %>% filter(probenID %in% sample_T1$probenID) %>% 
  left_join(ZirFlu$donorSamples) %>% select(-c(probenID, season, time)) %>%
  relocate(patientID)

View(flowCyto_data$`Baseline ex vivo`)


# check the association - linear --------------------
library(rstatix)

# proteins
protein_dat <- ZirFlu$protein_dat %>% select(selected_proteins$OlinkID) %>% 
  rownames_to_column("probenID") %>% filter(probenID %in% sample_T1$probenID) %>%
  left_join(ZirFlu$donorSamples %>% 
              filter(season == "2019") %>% select(probenID, patientID))

outcome <- list()
for (protein in selected_proteins$OlinkID) {
  outcome_temp <- c()
  
  for (cell in flowCyto_metadata$`Baseline ex vivo`$subcell_type) {
    dat_temp <- protein_dat %>% select(patientID, any_of(protein)) %>%
      full_join(flowCyto_data$`Baseline ex vivo` %>% select(c("Donor ID", any_of(cell))),
                by = c("patientID" = "Donor ID")) %>%
      column_to_rownames("patientID")
    
    cor_res <- cor_test(dat_temp)
    outcome_temp <- rbind(outcome_temp, cor_res)
  }
  outcome[[protein]] <- outcome_temp
}

# need to check all cells - no Th1, Th17, B cell show up
outcome_pval <- outcome %>% lapply(function(x) x %>% filter(p < 0.05))

plotList_1 <- list()
plotList_2 <- list()
for (celltype in flowCyto_metadata$`Baseline ex vivo`$subcell_type) {
  for (protein in selected_proteins$OlinkID) {
    plot_dat <- protein_dat %>% select(patientID, any_of(protein)) %>%
      full_join(flowCyto_data$`Baseline ex vivo` %>% select(c("Donor ID", any_of(celltype))),
                by = c("patientID" = "Donor ID")) %>%
      left_join(ZirFlu$donorInfo %>% filter(season == "2019"))
    
    plotList_1[[celltype]][[protein]] <- plot_dat %>% 
      ggscatter(x = celltype, y = protein, 
                color = "disease", add = "reg.line") +
      stat_cor(aes(color = disease), method = "pearson") +
      xlab(paste0(celltype, " percentage")) + 
      ylab(paste0("log2_expression_", 
                  ZirFlu$protein_annot$Assay[which(ZirFlu$protein_annot$OlinkID == protein) ]))
    
    plotList_2[[celltype]][[protein]] <- plot_dat %>% 
      ggscatter(x = celltype, y = protein, add = "reg.line") +
      stat_cor(method = "pearson") +
      xlab(paste0(celltype, " percentage")) + 
      ylab(paste0("log2_expression_", 
                  ZirFlu$protein_annot$Assay[which(ZirFlu$protein_annot$OlinkID == protein) ]))
    
  }
}

pdf(file = "temp/20230206_proteinT1_associCellType_T1.pdf", width = 15, onefile = TRUE)
for (celltype in flowCyto_metadata$`Baseline ex vivo`$subcell_type) {
  plot(cowplot::plot_grid(
    
    cowplot::plot_grid(plotList_1[[celltype]]$OID20477, 
                       plotList_1[[celltype]]$OID20481,
                       plotList_1[[celltype]]$OID20563,
                       plotList_1[[celltype]]$OID20611,
                       plotList_1[[celltype]]$OID20631,
                       nrow = 1),
    cowplot::plot_grid(plotList_2[[celltype]]$OID20477, 
                       plotList_2[[celltype]]$OID20481,
                       plotList_2[[celltype]]$OID20563,
                       plotList_2[[celltype]]$OID20611,
                       plotList_2[[celltype]]$OID20631,
                       nrow = 1),
    nrow = 2
  ))
}
dev.off()

plot_dat <- protein_dat %>% select(patientID, any_of(protein)) %>%
  full_join(flowCyto_data$`Baseline ex vivo` %>% select(c("Donor ID", any_of(celltype))),
            by = c("patientID" = "Donor ID")) %>%
  left_join(ZirFlu$donorInfo %>% filter(season == "2019"))

plot_dat %>% 
  ggscatter(x = celltype, y = protein, 
            color = "disease", add = "reg.line") +
  stat_cor(aes(color = disease), method = "pearson") +
  xlab(paste0(celltype, " percentage")) + 
  ylab(paste0("log2_expression_", 
              ZirFlu$protein_annot$Assay[which(ZirFlu$protein_annot$OlinkID == protein) ]))

plot_dat %>% 
  ggscatter(x = celltype, y = protein, add = "reg.line") +
  stat_cor(method = "pearson") +
  xlab(paste0(celltype, " percentage")) + 
  ylab(paste0("log2_expression_", 
              ZirFlu$protein_annot$Assay[which(ZirFlu$protein_annot$OlinkID == protein) ]))
