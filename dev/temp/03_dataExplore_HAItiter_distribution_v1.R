a<- ZirFlu$HAItiter_Feb2022 %>%
  select(condition, patientID, category, disease, 
         H1N1_titerFC:Byamagata.Phuket_titerFC) %>%
  pivot_longer(cols = H1N1_titerFC:Byamagata.Phuket_titerFC,
               names_to = "type",
               values_to = "log2_abTiter")
a %>% ggplot(aes(x = log2_abTiter, fill = type)) +
  geom_density(alpha = .3) +
  geom_vline(xintercept = 2, colour = "red", linetype = "longdash")
