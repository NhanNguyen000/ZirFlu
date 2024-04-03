ZirFlu$log2_cytokine <- ZirFlu$cytokine %>% replace(. == 0.01, 0) %>%
  get.log2() %>%
  rename(Condition = X...Condition) %>%
  mutate(Condition = ifelse(Condition == "Healthy", "Healthy", "Cirrhosis"))

ZirFlu$log2_cytokine %>% ggplot(aes(Measure, GRO....CXCL1., fill = Condition)) + 
  geom_boxplot() +
  scale_y_continuous(limits = c(0, 10))

ZirFlu$log2_cytokine %>% ggplot(aes(Measure, GRO....CXCL1., fill = Condition)) + 
  geom_violin() +
  geom_jitter(shape = 16, position = position_jitter(.02))

ggplot(ZirFlu$log2_cytokine, aes(Condition, GRO....CXCL1.)) + 
  geom_violin() +
  geom_jitter(shape = 16, position = position_jitter(0.02))

# boxplot link the data
cowplot::plot_grid(
  get.cytokin_violinPlot(ZirFlu$log2_cytokine, "Measure", "GRO....CXCL1."),
  get.cytokin_violinPlot(ZirFlu$log2_cytokine, "Measure", "IFN.."),
  get.cytokin_violinPlot(ZirFlu$log2_cytokine, "Measure", "IFN...1"),
  get.cytokin_violinPlot(ZirFlu$log2_cytokine, "Measure", "IL.6"),
  get.cytokin_violinPlot(ZirFlu$log2_cytokine, "Measure", "IL.8..CXCL8."),
  get.cytokin_violinPlot(ZirFlu$log2_cytokine, "Measure", "IL.9"),
  get.cytokin_violinPlot(ZirFlu$log2_cytokine, "Measure", "IL.10"),
  get.cytokin_violinPlot(ZirFlu$log2_cytokine, "Measure", "IL.13"),
  get.cytokin_violinPlot(ZirFlu$log2_cytokine, "Measure", "IL.18"),
  get.cytokin_violinPlot(ZirFlu$log2_cytokine, "Measure", "IL.21"),
  get.cytokin_violinPlot(ZirFlu$log2_cytokine, "Measure", "IL.22"),
  get.cytokin_violinPlot(ZirFlu$log2_cytokine, "Measure", "IP10..CXCL10."),
  get.cytokin_violinPlot(ZirFlu$log2_cytokine, "Measure", "RANTES..CCL5."),
  get.cytokin_violinPlot(ZirFlu$log2_cytokine, "Measure", "SDF.1.."),
  get.cytokin_violinPlot(ZirFlu$log2_cytokine, "Measure", "TNF.."),
  nrow = 3
)

# put everything in 1 plot
ggplot(ZirFlu$log2_cytokine, 
       aes(x = Condition, y = GRO....CXCL1., fill = Measure)) + 
  geom_violin() +
  geom_point(shape = 21, 
             position = position_jitterdodge(jitter.width = 0))
