rm(list = ls())

library(readxl)
library(tidyverse)

# functions used in the analysis -----------------------------------------------
get.log2 <- function(dat) {
  library(tidyverse)
  library(rstatix)
  outcome <- dat %>% mutate(across(where(is.numeric), ~log2(.x))) %>%
    mutate(across(where(is.numeric), ~ifelse(is.infinite(.x), 0, .x)))
  return(outcome)
}

# set up data list object and paths --------------------------------------------
ZirFlu <- list()

path_name <- "/vol/projects/BIIM/Influenza/ZirrFlu/metadata/"

# HAI antibody titers, use files from 11.2020 -----------------------------------

## HAI titer season 2019 -----------------------------------------------------

# HAI titer measurement
HAI_2019 <- read_excel(
  paste0(path_name, "20221101_Final_Info_HAI_Titer_Valerie&Nhan.xlsx"),
                          sheet = "2019-2020 HAI Titer") %>% 
  select(-c(15:19)) %>%
  drop_na("HAI Titer against Influenza A/H1N1") %>% as.data.frame()

names(HAI_2019) <- c("condition", "patientID", 
                     paste0(rep(c("H1N1_", "H3N2_", "Bvictoria_", "Byamagata_"), each = 3), 
                            HAI_2019[1,][3:14]))

# Foldchange of antibody (HAI titer) has been calculated
abFC_2019 <- read_excel(
  paste0(path_name, "20221101_Final_Info_HAI_Titer_Valerie&Nhan.xlsx"),
                            sheet = "2019-2020 Fold-increase") %>% 
  drop_na("Fold-increase") %>% as.data.frame()

names(abFC_2019) <- c("condition", "patientID", 
                            "H1N1_abFC", "H3N2_abFC", "Bvictoria_abFC", "Byamagata_abFC",
                            "vaccine_response")

# group the HAI measure and abFC
HAItiter_2019 <- HAI_2019 %>% slice(-1) %>% fill(condition) %>%
  full_join(abFC_2019 %>% 
              slice(-1) %>% fill(condition)) %>% 
  mutate(category = ifelse(vaccine_response == "non", "NR",
                           ifelse(vaccine_response %in% c("single", "double"), "Other", "TR"))) %>%
  relocate(condition, .after = last_col()) %>%
  rename_with(~sub("d0", "Baseline", .x)) %>% # Day 0 is baseline
  rename_with(~sub("d21-35", "Visit2", .x)) %>% # Day 21-35 is Visit 2
  rename_with(~sub("d55-75", "Visit3", .x)) %>% # Day 55 - 75 is Visit 3
  mutate_at(c(2:17), as.numeric) %>% get.log2()

rm(abFC_2019, HAI_2019)

## HAI titer season 2020 -----------------------------------------------------

# HAI titer measurement
HAI_2020 <- read_excel(paste0(path_name, "20221101_Final_Info_HAI_Titer_Valerie&Nhan.xlsx"),
                          sheet = "2020-2021 HAI Titer") %>% 
  select(-c(15:17)) %>%
  drop_na("HAI Titer against Influenza A/H1N1") %>% as.data.frame()

names(HAI_2020) <- c("condition", "patientID", 
                          paste0(rep(c("H1N1_", "H3N2_", "Bvictoria_", "Byamagata_"), each = 3),
                                 HAI_2020[1,][3:14]))

# Foldchange of antibody (HAI titer) has been calculated
abFC_2020 <- read_excel(paste0(path_name, "20221101_Final_Info_HAI_Titer_Valerie&Nhan.xlsx"),
                            sheet = "2020-2021 Fold-increase") %>% 
  drop_na("Fold-increase") %>% as.data.frame()

names(abFC_2020) <- c("condition", "patientID", 
                            "H1N1_abFC", "H3N2_abFC", "Bvictoria_abFC", "Byamagata_abFC",
                            "vaccine_response")

# group the HAI measure and abFC
HAItiter_2020 <- HAI_2020 %>% slice(-1) %>% fill(condition) %>%
  full_join(abFC_2020 %>% 
              slice(-1) %>% fill(condition)) %>%
  mutate(category = ifelse(vaccine_response == "non", "NR",
                           ifelse(vaccine_response %in% c("single", "double"), "Other", "TR"))) %>%
  relocate(condition, .after = last_col()) %>%
  rename_with(~sub("d0", "Baseline", .x)) %>% # Day 0 is baseline
  rename_with(~sub("d21-35", "Visit2", .x)) %>% # Day 21-35 is Visit 2
  rename_with(~sub("d56-104", "Visit3", .x)) %>% # Day 56-14 is Visit 3
  mutate_at(c(2:17), as.numeric) %>% get.log2() %>%
  filter(patientID != "Z-01-99-069")
# The participant Z-01-99-069 in season 2020 has azathioprine in his medication 
# and should be excluded from the analysis (info from the doctor).

rm(abFC_2020, HAI_2020)

# save data -------------------------------------------------------------------

ZirFlu$HAItiter <- HAItiter_2019 %>% mutate(season = "2019") %>%
  full_join(HAItiter_2020 %>% mutate(season = "2020")) %>%
  relocate("season", .after = "patientID") %>%
  mutate(condition = factor(condition, 
                            levels = c("Decompensated cirrhosis", 
                                       "Compensated cirrhosis", "Healthy"))) %>%
  mutate(disease = ifelse(condition == "Healthy", "Healthy", "Cirrhosis")) %>%
  mutate(disease = factor(disease,
                          levels = c("Cirrhosis", "Healthy")))


save(ZirFlu, file = "data/ZirFlu.RData")

