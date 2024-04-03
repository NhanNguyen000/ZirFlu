library(readxl)
library(tidyverse)
setwd("/vol/projects/BIIM/Influenza/ZirrFlu/metadata")

# HAI titer season 2019 -----------------------------------------------------
rawHAI_2019 <- read_excel("20221101_Final_Info_HAI_Titer_Valerie&Nhan.xlsx",
                      sheet = "2019-2020 HAI Titer")

cleanHAI_2019 <- rawHAI_2019 %>% select(-c(15:19)) %>%
  drop_na("HAI Titer against Influenza A/H1N1") %>% as.data.frame()

names(cleanHAI_2019) <- c("condition", "patientID", 
                          paste0(rep(c("H1N1_", "H3N2_", "Bvictoria_", "Byamagata_"), each = 3),
                                 cleanHAI_2019[1,][3:14]))


raw_abFC_2019 <- read_excel("20221101_Final_Info_HAI_Titer_Valerie&Nhan.xlsx",
                           sheet = "2019-2020 Fold-increase")

clean_abFC_2019 <- raw_abFC_2019 %>% 
  drop_na("Fold-increase") %>% as.data.frame()

names(clean_abFC_2019) <- c("condition", "patientID", 
                            "H1N1_abFC", "H3N2_abFC", "Bvictoria_abFC", "Byamagata_abFC",
                            "vaccine_response")


HAItiter_2019 <- cleanHAI_2019 %>% slice(-1) %>% fill(condition) %>%
  full_join(clean_abFC_2019 %>% slice(-1) %>% fill(condition))

# HAI titer season 2020 -----------------------------------------------------
rawHAI_2020 <- read_excel("20221101_Final_Info_HAI_Titer_Valerie&Nhan.xlsx",
                          sheet = "2020-2021 HAI Titer")

cleanHAI_2020 <- rawHAI_2020 %>% select(-c(15:17)) %>%
  drop_na("HAI Titer against Influenza A/H1N1") %>% as.data.frame()

names(cleanHAI_2020) <- c("condition", "patientID", 
                          paste0(rep(c("H1N1_", "H3N2_", "Bvictoria_", "Byamagata_"), each = 3),
                                 cleanHAI_2020[1,][3:14]))


raw_abFC_2020 <- read_excel("20221101_Final_Info_HAI_Titer_Valerie&Nhan.xlsx",
                            sheet = "2020-2021 Fold-increase")

clean_abFC_2020 <- raw_abFC_2020 %>% 
  drop_na("Fold-increase") %>% as.data.frame()

names(clean_abFC_2020) <- c("condition", "patientID", 
                            "H1N1_abFC", "H3N2_abFC", "Bvictoria_abFC", "Byamagata_abFC",
                            "vaccine_response")


HAItiter_2020 <- cleanHAI_2020 %>% slice(-1) %>% fill(condition) %>%
  full_join(clean_abFC_2020 %>% slice(-1) %>% fill(condition))

# remove unneeded variables -----------------------------------------------------
rm(raw_abFC_2019, rawHAI_2019, cleanHAI_2019, clean_abFC_2019,
   raw_abFC_2020, rawHAI_2020, cleanHAI_2020, clean_abFC_2020)


