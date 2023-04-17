library(tidyverse)
# cutoff -------
# (baseline < 20 ) == T1_low. (abFC>4) == abFC_high
# (baseline > 20) == T1_high. (abFC < 4) == abFC_low
baseline = log2(20)
log2_abFC = 4

# reclassification --------------
new_HAIgroups <- list()
new_HAIgroups$H1N1 <- ZirFlu$HAItiter %>%
  select(patientID, season, matches("H1N1"), vaccine_response, category, condition) %>%
  mutate(reclassify = ifelse(H1N1_T1 >= baseline, "H", "L")) %>%
  mutate(reclassify = paste0(reclassify, ifelse(H1N1_abFC >= log2_abFC, "H", "L")))

new_HAIgroups$H3N2 <- ZirFlu$HAItiter %>%
  select(patientID, season, matches("H3N2"), vaccine_response, category, condition) %>%
  mutate(reclassify = ifelse(H3N2_T1 >= baseline, "H", "L")) %>%
  mutate(reclassify = paste0(reclassify, ifelse(H3N2_abFC >= log2_abFC, "H", "L")))

new_HAIgroups$Bvictoria <- ZirFlu$HAItiter %>%
  select(patientID, season, matches("Bvictoria"), vaccine_response, category, condition) %>%
  mutate(reclassify = ifelse(Bvictoria_T1 >= baseline, "H", "L")) %>%
  mutate(reclassify = paste0(reclassify, ifelse(Bvictoria_abFC >= log2_abFC, "H", "L")))

new_HAIgroups$Byamagata <- ZirFlu$HAItiter %>%
  select(patientID, season, matches("Byamagata"), vaccine_response, category, condition) %>%
  mutate(reclassify = ifelse(Byamagata_T1 >= baseline, "H", "L")) %>%
  mutate(reclassify = paste0(reclassify, ifelse(Byamagata_abFC >= log2_abFC, "H", "L")))

# check groups
new_HAIgroups$H1N1 %>% count(season, reclassify)
new_HAIgroups$H3N2 %>% count(season, reclassify)
new_HAIgroups$Bvictoria %>% count(season, reclassify)
new_HAIgroups$Byamagata %>% count(season, reclassify)

## box plot --------------

strain <- "H1N1"
strain <- "H3N2"
strain <- "Bvictoria"
strain <- "Byamagata"
a <- new_HAIgroups[[strain]] %>% 
  select(patientID, season, matches(paste0(strain, "_T")), reclassify) %>% 
  gather(matches(strain), key = time, value = HAItiter) %>%
  arrange(season, patientID)
b <- a %>% filter(season == "2019")
b %>%
  ggplot(aes(time, HAItiter)) + geom_boxplot() + geom_point() + 
  geom_line(aes(group = patientID)) + theme_bw() +
  facet_wrap(vars(reclassify), nrow = 2)

b2 <- a %>% filter(season == "2020")
b2 %>%
  ggplot(aes(time, HAItiter)) + geom_boxplot() + geom_point() + 
  geom_line(aes(group = patientID)) + theme_bw() +
  facet_wrap(vars(reclassify), nrow = 2)
