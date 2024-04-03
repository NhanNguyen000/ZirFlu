ZirFlu$abTiters_Feb2022 %>% get.log2() %>%
  mutate(Condition = ifelse(Condition == "Healthy", "Healthy", "Cirrhosis")) %>%
  ggplot(aes(Condition, H1N1_d0)) + geom_violin() + 
  geom_jitter(shape = 16, position = position_jitter(0))

ZirFlu$abTiters_Feb2022 %>% get.log2() %>%
  mutate(Condition = ifelse(Condition == "Healthy", "Healthy", "Cirrhosis")) %>%
  ggplot(aes(Condition, H1N1_d21.35)) + geom_violin() + 
  geom_jitter(shape = 16, position = position_jitter(0))

# convert wide to long data
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
