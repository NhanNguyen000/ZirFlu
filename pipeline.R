## all the code is from NhanNguyen/ZirFLu_Rcodes.R
try(dev.off())
rm(list = ls())
source("R_functions_v1.R")

# 1. & 2. Load data into one ZirFlu list =======================================

# load the raw data
# use the lastest version of script stat with "01_loadData"
# excecute: run in the terminal "Rscript temp/01_loadData_[version].R" or
# setwd("script"); source("01_loadData_[version].R")
# outcome: rawData.Rdata 

# clean the raw data
load("rawDat.RData")
# use the lastest version of script stat with "02_cleanData"
# excecute: run in the terminal "Rscript temp/02_cleanData_[version].R" or 
# setwd("temp"); source("02_cleanData_[version].R")
# outcome: ZirFlu.Rdata 

# impute the protein data (if needed for WGCNA analysis)
load("ZirFlu.RData")
# script: 02_imputeData_v1.R --> impute missing value in protein data


# script: 02_filterMetabolite_v1.R

# 3. Explore data ==============================================================
# script: 03_dataExplore_cohortDescription_v1.R --> age, sex, and disease conditions
# script: 03_dataExplore_cohortDescription_v2.R --> Fig1.A (run in the terminal)

# script: 03_dataExplore_HAItiter_scatterPlot_v1.R --> scatter plot for antibody titers
# script: 03_dataExplore_HAItiter_boxplot_v1.R --> boxplot for antibody titers (with/out statistic test) --> Fig1.B
# script: 03_dataExplore_HAItiter_violinplot_v1.R

# script: 03_dataExplore_cytokine_v1.R

# script: 03_dataExplore_PCA_v1.R --> PCA plots for protein and metabolite -> Fig1.C
# script: 03_dataExplore_clusterTree_v1.R --> Cluster samples

# 4. Analyze data ==============================================================
## 4.1. T-test between healthy and cirhosiss group -----------------------------
# script: 04_ttest_allVariables_v1.R --> t-test for all HAI titer, cytokine, protein. metabolie to find candidate to look at
# script: 04_ttest_specificVariables_v1.R --> t-test and boxplot for specific candidate

## 4.2. lm() model for disease correction -----------------------------
# test with gender, age, condition correction: 04_lmModelCorrection --> 05_check_proteinCandidate_lmModel_v1.R

## 4.3. From padj metaboltie --> pathway --> function -----------------------------


## 4.4. WGCNA ------------------------------------------------------------------
# script: 04_wgcna_allSamples_v1.R
# script: 04_wgcna_perTime_v1.R -> Protein data - for T1 

# 5. Check interested variables ================================================
# script: 05_check_cytokineLinkAbTiter_v1.R --> check cytokine IL-15 for Saumya
# script: 05_check_proteinLinkAbTiter_v1.R --> check cytokine IL-15 (using Olinks - OID20562, uniprot: P40933) for Saumya
# script: 05_check_highCytokine_v1.R --> check sample with high cytokine value

# 6. Check the compatibility between data types ================================
# script: 06_check_cytokineLinkProtein_v1.R

# 7. Warp up plots =============================================================
cowplot::plot_grid(Fig1.A, Fig1.B, Fig1.C, 
                   labels = "AUTO", nrow = 1,
                   label_x = 0.02)
