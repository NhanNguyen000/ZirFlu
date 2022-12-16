# Chekc the identical participant between 2 seasons
a <- full_join(HAItiter_2019, HAItiter_2020, by = "patientID") %>%
  select(patientID, matches("vaccine")) %>% drop_na()