library(ggpubr)

cowplot::plot_grid(
  ZirFlu$log2_abTiters_Aug2022 %>% 
    ggplot(aes(x = donor, y = H1N1_titerFC, fill = donor)) + 
    geom_boxplot() + scale_y_continuous(limits = c(0, 10)) +
    ylab("Log2 of H1N1_titerFC") + theme_bw() +
    geom_hline(yintercept=2, linetype = "dashed") +
    theme(axis.text.x = element_text(angle = 45, hjust=1)),
  ZirFlu$log2_abTiters_Aug2022 %>% 
    ggplot(aes(x = donor, y = H3N2_titerFC, fill = donor)) + 
    geom_boxplot() + scale_y_continuous(limits = c(0, 10)) +
    ylab("Log2 of H3N2_titerFC") + theme_bw() +
    geom_hline(yintercept=2, linetype = "dashed") +
    theme(axis.text.x = element_text(angle = 45, hjust=1)),
  ZirFlu$log2_abTiters_Aug2022 %>% 
    ggplot(aes(x = donor, y = Bvictoria.Maryland_titerFC, fill = donor)) + 
    geom_boxplot() + scale_y_continuous(limits = c(0, 10)) +
    ylab("Log2 of Bvictoria.Maryland_titerFC") + theme_bw() +
    geom_hline(yintercept=2, linetype = "dashed") +
    theme(axis.text.x = element_text(angle = 45, hjust=1)),
  nrow = 1
)

Fig1.B <- cowplot::plot_grid(
  ggboxplot(ZirFlu$log2_abTiters_Feb2022, x = "Condition", y = "H1N1_titerFC",
            color = "Condition", paletter = "jco", add = "jitter") +
    stat_compare_means(method = "t.test"),
  ggboxplot(ZirFlu$log2_abTiters_Feb2022, x = "Condition", y = "H3N2_titerFC",
            color = "Condition", paletter = "jco", add = "jitter") +
    stat_compare_means(method = "t.test"),
  ggboxplot(ZirFlu$log2_abTiters_Feb2022, x = "Condition", y = "Bvictoria.Maryland_titerFC",
            color = "Condition", paletter = "jco", add = "jitter") +
    stat_compare_means(method = "t.test"),
  ggboxplot(ZirFlu$log2_abTiters_Feb2022, x = "Condition", y = "Byamagata.Phuket_titerFC",
            color = "Condition", paletter = "jco", add = "jitter") +
    stat_compare_means(method = "t.test"),
  nrow = 2
)
plot(Fig1.B)