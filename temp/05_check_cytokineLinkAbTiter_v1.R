#cytokine IL-15 (using Olinks - OID20562, uniprot: P40933) to antibody titer ------------------
cytokine_abTiters <- ZirFlu$log2_cytokine %>% 
  left_join(ZirFlu$log2_abTiters_Feb2022)
ZirFlu$log2_cytokine %>% filter(Measure == "d0") %>% select(IL.15) %>% unique()
ZirFlu$log2_cytokine %>% filter(Measure == "d21-35") %>% select(IL.15) %>% unique()

cytokine_abTiters %>% filter(Measure == "d0") %>% 
  ggplot(aes(x = H1N1_titerFC, y = IL.15, color = Category)) + 
  geom_point(position = position_jitterdodge(0.1)) + 
  geom_smooth(aes(color = "all"), method = "lm", se = FALSE)

coef(lm(IL.15 ~ H1N1_titerFC, data = cytokine_abTiters))

cytokine_abTiters %>% filter(Measure == "d0") %>% 
  ggplot(aes(x = H1N1_titerFC, y = IL.15, color = Category)) + 
  geom_point(position = position_jitterdodge(0.1)) + 
  geom_abline(intercept = 0.47529463, slope = 0.07419452)

cytokine_abTiters %>% filter(Measure == "d21-35") %>% 
  ggplot(aes(x = H1N1_titerFC, y = IL.15, color = Category)) + 
  geom_point(position = position_jitterdodge(0.1))

cytokine_abTiters %>%
  ggplot(aes(x = H1N1_titerFC, y = IL.15, color = Category)) + 
  geom_point(position = position_jitterdodge(0.1)) + facet_wrap(~Measure)