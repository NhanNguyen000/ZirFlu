library("ggalluvial")

metadata_flowchart <- ZirFlu$metadata2 %>% 
  select(Patient.ID, season, sex, age, condition, category) %>% distinct() %>%
  mutate(age = ifelse(age > 65, "old (>65)", "young"),
         gender = ifelse(sex == "m", "Male", "Female")) %>% 
  group_by(season, gender, age, condition, category) %>%
  count() %>% rename(donor = n, disease = condition) %>% 
  mutate(condition = ifelse(disease == "healthy", "healthy", "cirrhosis")) %>%
  mutate(age = factor(age, levels = c("young", "old (>65)")),
         disease = factor(disease, levels = c("healthy", "decompensated cirrhosis", 
                                              "compensated cirrhosis")),
         condition = factor(condition, levels = c("healthy", "cirrhosis"))) %>%
  relocate(donor, .after = last_col())


Fig1.A <- ggplot(metadata_flowchart, 
                 aes(y = donor, axis1 = gender, axis2 = age,
                     axis3 = disease, axis4 = condition)) +
  geom_alluvium(aes(fill = category), width = 8/12) +
  geom_stratum(width =8/12, alpha = 0.4, fill = "black") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("gender", "age", "disease", "condition"), expand = c(.10, .10)) +
  scale_fill_manual(values = c("dodgerblue", "darkolivegreen2", "indianred2")) +
  theme_bw() + theme(legend.position = "top")
plot(Fig1.A)
