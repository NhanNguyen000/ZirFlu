ZirFlu$HAItiter_Feb2022 %>% 
  ggplot(aes(condition, H1N1_d0)) + geom_violin() + 
  geom_jitter(shape = 16, position = position_jitter(0))

a<- ZirFlu$HAItiter_Feb2022 %>%
  select(condition, patientID, category, disease, 
         H1N1_titerFC:Byamagata.Phuket_titerFC) %>%
  pivot_longer(cols = H1N1_titerFC:Byamagata.Phuket_titerFC,
               names_to = "type",
               values_to = "log2_abTiter")
a %>% ggplot(aes(x = type, y = log2_abTiter, fill = type)) + 
  geom_violin() +
  geom_point(shape = 16, 
             position = position_jitterdodge(jitter.height = 0.2, jitter.width = 0.2)) 

a %>% ggplot(aes(x = condition, y = log2_abTiter, fill = condition)) + 
  geom_violin() +
  geom_point(shape = 16, 
             position = position_jitterdodge(jitter.height = 0.2, jitter.width = 0.2)) +
  facet_wrap(~type) +
  geom_hline(yintercept = 2, colour = "red", linetype = "longdash")


# convert wide to long data ================ old codes
a <- ZirFlu$abTiters_Feb2022 %>% 
  mutate(Condition = ifelse(Condition == "Healthy", "Healthy", "Cirrhosis")) %>%
  pivot_longer(cols = 3:10, names_to = "ab", values_to = "log2ab_level") %>%
  separate(ab, c("ab_type", "time"), "_") %>%
  get.log2()

a %>% ggplot(aes(x = time, y = log2ab_level, fill = Condition)) + 
  geom_violin() +
  # geom_point(shape = 16, position = position_jitterdodge(jitter.width = 0)) +
  facet_wrap(~ab_type)


a %>% ggplot(aes(x = interaction(time, Condition), y = log2ab_level, fill = Condition)) + 
  geom_violin() +
  geom_line(aes(group = interaction(PatientID, Condition))) +
  facet_wrap(ab_type~.)
