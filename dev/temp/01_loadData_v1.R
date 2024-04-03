library(OlinkAnalyze)
library(openxlsx)

setwd("/vol/projects/CIIM/Influenza/ZirrFlu/NhanNguyen")
rawDat <- list()

rawDat$metadata <- read.table(file = '../metadata/zirrflu_meta.csv', sep = '\t', header=T)
rawDat$HAItiter_Feb2022 <- read.table("../metadata/HAI-Titer_Data_from202202_NN.csv", 
                                      sep = ";", header = TRUE)
rawDat$HAItiter_Aug2022 <- read.table("../metadata/HAI-Titer_Data_from202208_NN.csv", 
                                      sep = ";", header = TRUE)
rawDat$cytokine_Feb2022 <- read.table("../metadata/Cytokine_Data_from202202_NN.csv", 
                              sep = ";", header = TRUE, dec = ",")
rawDat$protein <- read_NPX(filename = "../proteomic/raw_data/20212645_Li_NPX_2022-02-02.csv")
rawDat$metabolite <- read.xlsx(xlsxFile = '../metabolic/raw_data/spreadsheets/DATA.xlsx',
                               sheet = 'ions')
rawDat$metabolite_annot <- read.xlsx(xlsxFile = '../metabolic/raw_data/spreadsheets/DATA.xlsx',
                               sheet = 'annotation')
setwd("/vol/projects/CIIM/Influenza/ZirrFlu/ZirFlu_NhanNguyen")
save(rawDat, file = "rawDat.RData")
