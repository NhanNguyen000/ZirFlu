## age, sex, and disease conditions ----------------
df <- ZirFlu$metadata %>% 
  select(Patient.ID, sex, age, HEP_HLA_B27, medication, Disease, child_pugh_score, 
         cirrhosis.ethiology, cirrhosis.therapy, condition) %>%
  remove_rownames() %>% distinct() %>%
  mutate_if(is.character, as.factor)

ggplot(df, aes(x=age, group = sex, fill = sex)) + 
  geom_histogram(binwidth = 2) + facet_wrap(~sex)
ggplot(df) + geom_bar(aes(x=condition, fill = sex))