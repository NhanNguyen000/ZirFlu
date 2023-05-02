# Metadata: clean, check, calculate  ==============================
get.vaccine_resp <- function(dat, pre_vac, post_vac) {
  outcome <- as.data.frame(dat$Patient.ID)
  
  for (i in 1:nrow(dat)) {
    outcome[i, paste0(pre_vac, "_FC")] <- ifelse(dat[i, pre_vac] >=10, dat[i, post_vac]/dat[i, pre_vac],
                                   ifelse(dat[i, post_vac] > 40, dat[i, post_vac], NA))
    outcome[i, paste0(pre_vac, "_vacResponse")] <- ifelse(outcome[i, paste0(pre_vac, "_FC")] >= 4, TRUE, FALSE)
  }
  return(outcome)
}

get.log2 <- function(dat) {
  library(tidyverse)
  library(rstatix)
  outcome <- dat %>% mutate(across(where(is.numeric), ~log2(.x))) %>%
    mutate(across(where(is.numeric), ~ifelse(is.infinite(.x), 0, .x)))
  return(outcome)
}

get.completeCase_dat <- function(dat, time_point) {
  library(tidyverse)
  patients <- intersect(dat[which(dat$Time == time_point[1]),]$Patient.ID,
                        dat[which(dat$Time == time_point[2]),]$Patient.ID)
  outcome <- dat %>% filter(Patient.ID %in% patients) %>%
    filter(Time %in% time_point)
  return(outcome)
}

check.NAs <- function(dat) {
  library(naniar)
  library(visdat)
  library(tidyverse)
  
  print(paste0("The number of missing values: ", n_miss(dat)))
  print(paste0("The number of complete values: ", n_complete(dat)))
  
  cat("The percents of missingness per variable - range: ")
  cat(dat %>% miss_var_summary() %>% select(pct_miss) %>% range())
  
  cat("\nThe percents of missingness per case/observation - range: ")
  cat(dat %>% miss_case_summary() %>% select(pct_miss) %>% range())
  
  vis_miss(dat, cluster = TRUE, sort_miss = TRUE)
  
}

test.normality <- function(dat) {
  test_temp <- list()
  for(idx in names(dat)) {
    test_temp[[idx]] <- shapiro.test(dat[, idx]) %>% unlist()
  }
  outcome <- as.data.frame(matrix(unlist(test_temp), ncol = 4, byrow = TRUE))
  names(outcome) <- names(test_temp[[idx]])
  rownames(outcome) <- names(test_temp)
  return(outcome[1:3])
}

# Plots -  statistic description ===============================================
get.pca_plot <- function(pca, metadat, groupType) {
  library(tidyverse)
  library(caret)
  
  ggplot(data.frame(pca$x), aes(PC1, PC2, color = metadat[, groupType])) +
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

get.outlier_inPCA <- function(pca, outliers) {
  library(tidyverse)
  library(caret)
  labels <- ifelse(rownames(pca$x) %in% outliers, rownames(pca$x), "")
  colored_point <- ifelse(labels != "", "red", "black")
  
  ggplot(data.frame(pca$x), aes(PC1, PC2, color = colored_point)) + 
    geom_point() + theme_classic() + scale_color_identity() +
    labs(x = paste0("PC1 [", round(100*summary(pca)$importance["Proportion of Variance", "PC1"], 
                                   digits = 2), "%]"),
         y = paste0("PC2 [", round(100*summary(pca)$importance["Proportion of Variance", "PC2"], 
                                   digits = 2), "%]")) + 
    geom_text(aes(label = labels), colour = "red", hjust = 0, vjust = 1.5)  
}

get.comparMeasures_plot <- function(log2_dat, x_name, y_name, color_group, abline_val) {
  library(ggplot2)
  
  ggplot(log2_dat, aes(x = log2_dat[, x_name], y = log2_dat[, y_name],
                       color = log2_dat[, color_group])) + 
    geom_point(size = 3,
               position = position_jitterdodge(jitter.width = .2, jitter.height = .2)) +
    xlab(paste0("Log2 of ", x_name)) + ylab(paste0("Log2 of ", y_name)) +
    scale_color_discrete(name = color_group) +
    xlim(-1, 15) + ylim(-1, 15) +
    geom_vline(xintercept = log2(abline_val)) +
    geom_hline(yintercept = log2(abline_val)) +
    geom_text(aes(x = 10, y = log2(abline_val)+0.5,
                  label = paste0("log2(", abline_val, ")"))) +
    geom_text(aes(x = log2(abline_val)-0.5 , y = 10,
                  label =paste0("log2(", abline_val, ")")),
              angle=90)
}

get.scatter_plot <- function(log2_dat, x_name, y_name, color_group, shape_group) {
  library(ggplot2)
  
  ggplot(log2_dat, aes(x = log2_dat[, x_name], y = log2_dat[, y_name], 
                 color = log2_dat[, color_group],
                 shape = log2_dat[, shape_group])) + 
    geom_point(size = 3, 
               position = position_jitterdodge(jitter.width = .2, jitter.height = .2)) + 
    xlab(paste0("Log2 of ", x_name)) + ylab(paste0("Log2 of ", y_name)) +
    scale_color_discrete(name = color_group) +
    scale_shape_discrete(name = shape_group) +
    xlim(-1, 12) + ylim(0, 15)
}

get.scatterPlot2 <- function(log2_dat, x_name, y_name, color_group) {
  library(ggplot2)
  
  ggplot(log2_dat, aes(x = log2_dat[, x_name], y = log2_dat[, y_name], 
                       color = log2_dat[, color_group])) + 
    geom_point(position = position_jitterdodge(jitter.width = .2, jitter.height = .2)) + 
    xlab(paste0("Log2 of ", x_name)) + ylab(paste0("Log2 of ", y_name)) +
    scale_color_discrete(name = color_group)
}

get.scatterPlot3 <- function(log2_dat, x_name, y_name, color_group) {
  library(ggplot2)
  
  ggplot(log2_dat, aes(x = log2_dat[, x_name], y = log2_dat[, y_name])) + 
    geom_point(aes(color = log2_dat[, color_group]),
               position = position_jitterdodge(0.1)) +
    geom_smooth(method = "lm", se = FALSE) +
    xlab(paste0("Log2 of ", x_name)) + 
    ylab(paste0("Log2 of ", y_name)) +
    scale_color_manual(values = c("TR" = "indianred2",
                                  "Other" = "darkolivegreen2",
                                  "NR" = "dodgerblue"),
                       name = color_group) +
    theme_bw()
}
get.scatterPlot4 <- function(log2_dat, x_name, y_name, color_group) {
  library(ggplot2)
  
  ggplot(log2_dat, aes(x = log2_dat[, x_name], y = log2_dat[, y_name])) + 
    geom_point(aes(color = log2_dat[, color_group]),
               position = position_jitterdodge(0.1)) +
    geom_smooth(method = "lm", se = FALSE) +
    xlab(paste0("Log2 of ", x_name)) + 
    ylab(paste0("Log2 of ", y_name)) + labs(color = color_group) +
    theme_bw()
}


get.cytokin_violinPlot <- function(log2_dat, x_name, y_name) {
  library(ggplot2)
  
  log2_dat %>% 
    ggplot(aes(x = log2_dat[, x_name], y = log2_dat[, y_name])) + 
    geom_violin() + geom_jitter(shape = 16, position = position_jitter(0.02)) + 
    xlab(paste0(x_name)) + ylab(paste0("Log2 of ", y_name)) +
    geom_line(aes(group = PatientID)) + facet_wrap(~Condition)
}

get.boxplot <- function(dat, y_value, ylim_value = NA) {
  library(ggplot2)
  ggplot(dat, aes(x = dat[, "category"], y = dat[, y_value], 
                  fill = dat[, "condition"])) +
    geom_boxplot() + ylim(0, ylim_value) +
    geom_jitter(color="black", size=0.5, alpha=0.9)
}

# Statistic tests ==============================================================
get.ttest <- function(dat, selected_var) {
  library(tidyverse)
  library(rstatix)
  outcome <- as.data.frame(matrix(nrow = 0, ncol = 8))
  names(outcome) <- c(".y.", "group1", "group2", "n1", "n2","statistic", "df", "p")
  
  for (i in selected_var) {
    dat_temp <- dat %>% select(Patient.ID, group, {{i}}) %>% rename(value = {{i}})
    res_temp <- dat_temp %>% t_test(value ~ group)
    res_temp[,1] <- i
    outcome <- rbind(outcome, res_temp)
  }
  return(outcome)
}

get.pairedTtest <- function(dat, list_vars) {
  library(tidyverse)
  library(rstatix)
  outcome <- as.data.frame(matrix(nrow = 0, ncol = 9))
  names(outcome) <- c("condition", ".y.", "group1", "group2", "n1", "n2","statistic", "df", "p")
  for (i in list_vars) {
    dat_temp <- dat %>% 
      select(Patient.ID, group, condition, {{i}}) %>% rename(value = {{i}})
    res_temp <- dat_temp %>% group_by(condition) %>% 
      t_test(value ~ group, paired = TRUE)
    res_temp[,2] <- i
    
    outcome <- rbind(outcome, res_temp)
  }
  return(outcome)
}

get.lmTest_correction <- function(variableList, abTiter, inputDat) {
  outcome <- data.frame(targetVariable = character(), 
                        independentVariable = character(),
                        estimate = double(), std.error = double(),
                        statistic = double(), p.value = double())
  for (val in variableList) {
    dat_temp <- inputDat %>% 
      select(all_of(val), sex, age, all_of(abTiter), disease) %>%
      rename(Variable = {{val}}, abTiter = {{abTiter}})
    
    lm_result <- tidy(lm(Variable ~sex + age + abTiter + disease, data = dat_temp)) %>%
      select(term, estimate, std.error, statistic, p.value) %>%
      as.data.frame()
    outcome[nrow(outcome) + 1, ] <- c(val, lm_result[4, ])
    outcome[nrow(outcome) + 1, ] <- c(val, lm_result[5, ])
  }
  return(outcome)
}

get.lmTest_correction_v2 <- function(abFC, variableList, inputDat) {
  outcome <- data.frame(targetVariable = character(), 
                        independentVariable = character(),
                        estimate = double(), std.error = double(),
                        statistic = double(), p.value = double())
  for (protein in variableList) {
    dat_temp <- inputDat %>% 
      select(all_of(abFC), sex, age, all_of(protein), disease) %>%
      rename(abFC = {{abFC}}, proteinVal = {{protein}})
    
    lm_result <- tidy(lm(abFC ~ sex + age + proteinVal + disease, data = dat_temp)) %>%
      select(term, estimate, std.error, statistic, p.value) %>%
      as.data.frame()
    
    outcome[nrow(outcome) + 1, ] <- c(protein, lm_result[4, ])
    outcome[nrow(outcome) + 1, ] <- c(protein, lm_result[5, ])
  }
  return(outcome)
}

get.var_onlyAbTiter <- function(var_associAbTiter, var_associDisease, abTiter) {
  vars <- setdiff(var_associAbTiter[[abTiter]]$targetVariable,
                      var_associDisease[[abTiter]]$targetVariable)
  return(vars)
}

# WGCNA =====================================================================
get.WGCNA_sampleTree <- function(dat) {
  sampleTrees <- hclust(dist(dat), method = "complete")
  plot(sampleTrees,
       main = "Sample clustering to detect outliers", 
       sub="", xlab="", cex.lab = 1.5,
       cex.axis = 1.5, cex.main = 2)
  return(sampleTrees)
}

rm.samples <- function(sampleTree, input_dat, cutoff) {
  plot(sampleTree,
       main = "Sample clustering to remove outliers", 
       sub="", xlab="", cex.lab = 1.5,
       cex.axis = 1.5, cex.main = 2)
  abline(h = cutoff, col = "red")
  
  clust <- cutreeStatic(sampleTree, cutHeight =  cutoff, minSize = 10)
  keepSamples <- (clust == 1) # clust 1 contains the samples to keep
  
  output_dat <- input_dat[keepSamples, ]
  return(output_dat)
}

get.WGCNA_powerTables <- function(dat) {
  powers = c(c(1:10), seq(12,20, by = 2))
  powerTable <- pickSoftThreshold(dat, powerVector = powers, verbose = 5)[[2]]
  
  plotCols = c(2,5,6,7)
  colNames = c("Scale Free Topology Model Fit", "Mean connectivity", "Median connectivity", "Max connectivity")
  par(mfcol = c(2,2))
  for (col in 1:length(plotCols)) {
    plot(powerTable[,1], -sign(powerTable[,3]) * powerTable[,2],
         xlab = "Soft threshold (power)", 
         ylab = colNames[col], ylim = c(0, max(powerTable[,plotCols[col]])) ,
         type = "n", main = colNames[col])
    addGrid()
    if (col == 1) text(powerTable[,1], -sign(powerTable[,3]) * powerTable[,2], labels = powers)
    else text(powerTable[,1], powerTable[,plotCols[col]], labels = powers)
  }
}


get.WGCNAnet <- function(input_dat, selected_power, netType) {
  net <- blockwiseModules(input_dat, power = selected_power, TOMType = netType, 
                          minModuleSize = 30, reassignThreshold = 0,
                          mergeCutHeight = 0.25, numericLabels = TRUE,
                          pamRespectsDendro = FALSE, saveTOMs = FALSE,
                          verbose = 3)
  mergedColors <- labels2colors(net$colors)
  plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                      "Module colors", dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05,
                      main = "Consensus dendrogram and module colors")
  return(net)
}

get.module_elements <- function(WGCNAnet) {
  moduleColors <- labels2colors(WGCNAnet$colors)
  outcome <- list()
  for (moduleCol in unique(moduleColors)) {
    outcome[[moduleCol]] <- names(WGCNAnet$colors[which(moduleColors == moduleCol)])
  }
  return(outcome)
}

get.intraModule_connectivity <- function(input_dat, WGCNAnet, WGCNA_modulesAnnot) {
  connectTable <- intramodularConnectivity.fromExpr(input_dat, colors = WGCNAnet$colors)
  rownames(connectTable) <- names(input_dat)
  
  connectModule <- list()
  for (name in names(WGCNA_modulesAnnot)) {
    connectModule[[name]] <- connectTable %>% rownames_to_column("OlinkID") %>%
      filter(OlinkID %in% WGCNA_modulesAnnot[[name]]$OlinkID) %>%
      full_join(WGCNA_modulesAnnot[[name]])
  }
  return(connectModule)
}

get.intraModule_connectivity_metabolites <- function(input_dat, WGCNAnet, WGCNA_modulesAnnot) {
  connectTable <- intramodularConnectivity.fromExpr(input_dat, colors = WGCNAnet$colors)
  rownames(connectTable) <- names(input_dat)
  
  connectModule <- list()
  for (name in names(WGCNA_modulesAnnot)) {
    connectModule[[name]] <- connectTable %>% rownames_to_column("ionIdx") %>%
      mutate(ionIdx = as.numeric(ionIdx)) %>%
      filter(ionIdx %in% WGCNA_modulesAnnot[[name]]$ionIdx) %>%
      full_join(WGCNA_modulesAnnot[[name]])
  }
  return(connectModule)
}

get.moduleTrait_relation <- function(MEs, datTraits) {
  modTraitCor <- cor(MEs, datTraits, use = "p")
  modTraitPvalue <- corPvalueStudent(modTraitCor, nrow(datTraits))
  
  textMatrix <- paste0(signif(modTraitCor, 2), "\n(", signif(modTraitPvalue, 2), ")")
  dim(textMatrix) <- dim(modTraitCor)
  
  sizeGrWindow(10,6); par(mar = c(6, 8.5, 3, 3))
  labeledHeatmap(modTraitCor, xLabels = names(datTraits), 
                 yLabels = names(MEs), ySymbols = names(MEs),
                 colorLabels = FALSE, colors = blueWhiteRed(50),
                 textMatrix = textMatrix, setStdMargins = FALSE,
                 cex.text = 1, zlim = c(-1, 1),
                 main = "Module - trait relationships")
  return(list("modTraitCor" = modTraitCor,
              "modTraitPvalue" = modTraitPvalue))
}

get.geneModuleMembership <- function(MEs, input_dat) {
  modNames <- substring(names(MEs), 3)
  
  geneMM <- as.data.frame(cor(input_dat, MEs, use = "p"))
  names(geneMM) <- paste0("MMcor", modNames)
  
  geneMM_pvalue <- as.data.frame(corPvalueStudent(as.matrix(geneMM), nrow(input_dat)))
  names(geneMM_pvalue) <- paste0("p.MMcor", modNames)
  return(list("geneMM" = geneMM, "geneMM_pvalue" = geneMM_pvalue))
}

get.geneTrait_relation <- function(input_dat, datTrait) {
  geneTrait <- as.data.frame(cor(input_dat, datTrait, use = "p"))
  geneTrait_pvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTrait), nrow(input_dat)))
  outcome <- cbind(geneTrait, geneTrait_pvalue)
  colnames(outcome) <- c("GScor", "p.GScor")
  return(outcome)
}

get.geneModuleTrait_relation <- function(module, WGCNA_network,
                                         geneModuleMembership, 
                                         geneTraitSignificance) {
  moduleGenes <- labels2colors(WGCNA_network$colors) == module
  
  verboseScatterplot(abs(geneModuleMembership[moduleGenes,
                                              grep(module, names(geneModuleMembership))]),
                     abs(geneTraitSignificance[moduleGenes, 1]),
                     xlab = paste0("Module Membership in ", module, " module"),
                     ylab = "Gene significance for selected trait",
                     main = paste("Module membership vs. gene significance \n"),
                     cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
  
}

get.geneModuleTrait_relation2 <- function(module, WGCNA_network,
                                         geneModuleMembership, 
                                         geneTraitSignificance) {
  moduleGenes <- labels2colors(WGCNA_network$colors) == module
  
  verboseScatterplot(geneModuleMembership[moduleGenes,
                                              grep(module, names(geneModuleMembership))],
                     geneTraitSignificance[moduleGenes, 1],
                     xlab = paste0("Gene - module Membership in ", module, " module"),
                     ylab = "Gene-trait correlation for selected trait",
                     main = paste("Gene - module membership vs. gene-trait correlation \n"),
                     cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
  
}

get.geneIntraConnect_vsTrait <- function(module,intraMod_connect, geneTraitCor) {
  verboseScatterplot(intraMod_connect[[module]]$kWithin,
                     geneTraitCor[intraMod_connect[[module]]$OlinkID, 1],
                     xlab = "Intra module connectivity",
                     ylab = "Gene-trait correlation for selected trait",
                     main = paste("Intra module connectivity vs. gene-trait correlation \n"),
                     cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
  
}

# Annotation =====================================================================
## protein annotation ----------------------------------------------------
get.proteinAnnot <- function(annot, input) {
  outcome <- annot %>% filter(OlinkID %in% input)
  return(outcome)
}

## metabolite annotation ----------------------------------------------------
get.keggBrite <- function(kegg_dat) {
  brite_temp <- kegg_dat$BRITE
  names(brite_temp) <- kegg_dat$ENTRY
  
  keggBrite <- list()
  for (id in names(brite_temp)) keggBrite[[id]] <- brite_temp[[id]]
  return(keggBrite)
}

get.keggTaxonomy <- function(keggBrite) {
  library(tidyverse)
  outcome <- tibble(V = keggBrite) %>% unnest_wider(V, names_sep = "") %>%
    mutate(id = names(keggBrite)) %>% column_to_rownames("id")
  return(outcome)
}

get.brite_class <- function(brite_case) {
  briteIds <- grep("BR:br", brite_case)
  brite_temp <- list()
  
  for (idx in 1: length(briteIds)) {
    br_name <- brite_case[briteIds[idx]]
    
    if (idx == length(briteIds)) {
      brite_temp[[br_name]] <- brite_case[briteIds[idx] : length(brite_case)]
    } else {
      brite_temp[[br_name]] <- brite_case[briteIds[idx] : (briteIds[idx+1]-1)]
    }
  }
  return(brite_temp)
}


get.kegg_briteGroups <- function(keggBrite) {
  
  keggBrite_subgroup <- list()
  for (kegg_id in names(keggBrite)) {
    keggBrite_subgroup[[kegg_id]] <- get.brite_class(keggBrite[[kegg_id]])
  }
  return(keggBrite_subgroup)
}

get.aBriteGroup <- function(kegg_briteGroups, groupName) {
  brite_group <- list()
  for(kegg_id in names(kegg_briteGroups)) {
    if (groupName %in% names(kegg_briteGroups[[kegg_id]])) {
      brite_group [[kegg_id]] <- kegg_briteGroups[[kegg_id]][[groupName]]
    }
  }
  return(brite_group)
}

get.briteNote <- function(keggTaxonomy) {
  library(tidyverse)
  library(stringr)
  
  if (ncol(keggTaxonomy) <= 5) outcome <- keggTaxonomy %>% mutate_if(is.character, str_trim)
  else {outcome <- keggTaxonomy %>% 
    mutate(V6 = ifelse(is.na(V6), NA, "Have multiple name, please check")) %>% 
    select(V1:V6) %>% mutate_if(is.character, str_trim)}
  
  return(outcome)
}

get.metaboliteAnnot <- function(annot, input) {
  outcome <- annot %>% filter(ionIdx %in% input)
  return(outcome)
}



