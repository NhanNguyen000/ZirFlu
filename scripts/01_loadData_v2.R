library(OlinkAnalyze)
library(openxlsx)

setwd("/vol/projects/CIIM/Influenza/ZirrFlu/NhanNguyen")
rawDat <- list()

## HAI antibody titers ---------------------------------------------------------
# run the script stated with: 00_cleanHAI_files2020Nov
rawDat$HAItiter_2019 <- HAItiter_2019
rawDat$HAItiter_2020 <- HAItiter_2020

## metadata --------------------------------------------------------------------
rawDat$metadata_2019 <- read.table(file = '../metadata/zirrflu_meta.csv', sep = '\t', header=T)

## load other data types -------------------------------------------------------
rawDat$cytokine_Feb2022 <- read.table("../metadata/Cytokine_Data_from202202_NN.csv", 
                                      sep = ";", header = TRUE, dec = ",")
rawDat$protein <- read_NPX(filename = "../proteomic/raw_data/20212645_Li_NPX_2022-02-02.csv")
rawDat$metabolite <- read.xlsx(xlsxFile = '../metabolic/raw_data/spreadsheets/DATA.xlsx',
                               sheet = 'ions')
rawDat$metabolite_annot <- read.xlsx(xlsxFile = '../metabolic/raw_data/spreadsheets/DATA.xlsx',
                                     sheet = 'annotation')

## save data -------------------------------------------------------------------
setwd("/vol/projects/CIIM/Influenza/ZirrFlu/ZirFlu_NhanNguyen")
save(rawDat, file = "rawDat.RData")
