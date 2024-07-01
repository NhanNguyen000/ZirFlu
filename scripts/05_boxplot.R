rm(list = ls())

library(tidyverse)
library(ggpubr)

# NOTE: The color code for each condition changed a bit, because other people saids the contrast is not clear in these specific plots 

# load the data list object------------------------------------------------------
load("data/ZirFlu.RData")

metadata <- ZirFlu$donorSamples %>% 
  full_join(ZirFlu$donorInfo) %>%
  drop_na() %>%
  mutate(condition = factor(condition, 
                            levels = c("Healthy", "Compensated cirrhosis", "Decompensated cirrhosis")))

# prepare the data -----------------------------------------------------------------
metaboliteDat <- metadata %>% 
  left_join(ZirFlu$metabolite_dat %>% 
              t() %>% as.data.frame() %>%
              rownames_to_column("ionIdx") %>% mutate(ionIdx = as.numeric(ionIdx)) %>%
              left_join(ZirFlu$metabolite_annot %>% select(ionIdx, Formula) %>% distinct()) %>%
              select(-ionIdx) %>% column_to_rownames("Formula") %>% 
              t() %>% as.data.frame() %>%
              rownames_to_column("probenID"))

metaboliteDat_plot <- metaboliteDat %>% filter(time %in% c("Baseline", "Visit2")) %>%
  full_join(ZirFlu$HAItiter)

# plotting  -----------------------------------------------------------------
my_comparisons <- list(c("Healthy", "Compensated cirrhosis"), 
                       c("Healthy", "Decompensated cirrhosis"), 
                       c("Compensated cirrhosis", "Decompensated cirrhosis"))

## boxplot --------------------------------
year <- "2019"
year <- "2020"

boxplot_mebo <- cowplot::plot_grid(
  metaboliteDat_plot %>% filter(season == year) %>%
    ggboxplot(x = "condition", y = "C20H34O5", 
              color = "condition", palette = c("#1515FF", "#CECE0E", "#918D8D"),
              add = "jitter", add.params = list(size = 3)) + rotate_x_text(25) +
    stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "t.test") +
    facet_wrap("time"),
  metaboliteDat_plot %>% filter(season == year) %>%
    ggboxplot(x = "condition", y = "C5H11NO2", 
              color = "condition", palette = c("#1515FF", "#CECE0E", "#918D8D"),
              add = "jitter", add.params = list(size = 3)) + rotate_x_text(25) +
    stat_compare_means(comparisons = my_comparisons, method = "t.test", label = "p.signif") +
    facet_wrap("time"),
  nrow = 1
)

# Save the plot as SVG
boxplot_mebo  
ggsave(paste0("output/boxplot_mebo _season", year, ".svg"), 
       boxplot_mebo , device = "svg")
ggsave(paste0("output/boxplot_mebo_season", year, ".png"), 
       boxplot_mebo , device = "png")

## scatter plot season 2019 --------------------------------
year <- "2019"
year <- "2020"

plotDat <- metaboliteDat_plot %>% filter(season == year) %>%
  rename("H1N1" = "H1N1_Visit2", "H3N2" = "H3N2_Visit2",
         "B/Victoria" = "Bvictoria_Visit2", "B/Yamagata" = "Byamagata_Visit2")

# for "C20H34O5"
C20H34O5_scaterplot <- cowplot::plot_grid(
  plotDat %>%
    ggscatter(x = "H1N1", y = "C20H34O5", 
              color = "condition", palette = c("#1515FF", "#CECE0E", "#918D8D"),
              add = "reg.line", xlim = c(5, 14)) +
    stat_cor(aes(color = condition), method = "spearman", label.x = 5),
  plotDat %>%
    ggscatter(x = "H3N2", y = "C20H34O5",
              color = "condition", palette = c("#1515FF", "#CECE0E", "#918D8D"),
              add = "reg.line", xlim = c(5, 14)) + 
    stat_cor(aes(color = condition), method = "spearman", label.x = 5),
  plotDat %>%
    ggscatter(x = "B/Victoria", y ="C20H34O5", 
              color = "condition", palette = c("#1515FF", "#CECE0E", "#918D8D"),
              add = "reg.line", xlim = c(5, 14)) +
    stat_cor(aes(color = condition), method = "spearman", label.x = 5),
  plotDat %>%
    ggscatter(x = "B/Yamagata", y = "C20H34O5", 
              color = "condition", palette = c("#1515FF", "#CECE0E", "#918D8D"),
              add = "reg.line", xlim = c(5, 14)) +
    stat_cor(aes(color = condition), method = "spearman", label.x = 5),
  nrow = 1
)

# for "C5H11NO2"
C5H11NO2_scaterplot <- cowplot::plot_grid(
  plotDat %>%
    ggscatter(x = "H1N1", y = "C5H11NO2", 
              color = "condition", palette = c("#1515FF", "#CECE0E", "#918D8D"),
              add = "reg.line", xlim = c(5, 14)) +
    stat_cor(aes(color = condition), method = "spearman", label.x = 8),
  plotDat %>%
    ggscatter(x = "H3N2", y = "C5H11NO2", 
              color = "condition", palette = c("#1515FF", "#CECE0E", "#918D8D"),
              add = "reg.line", xlim = c(5, 14)) +
    stat_cor(aes(color = condition), method = "spearman", label.x = 8),
  plotDat %>%
    ggscatter(x = "B/Victoria", y ="C5H11NO2", 
              color = "condition", palette = c("#1515FF", "#CECE0E", "#918D8D"),
              add = "reg.line", xlim = c(5, 14)) +
    stat_cor(aes(color = condition), method = "spearman", label.x = 8),
  plotDat %>%
    ggscatter(x = "B/Yamagata", y = "C5H11NO2", 
              color = "condition", palette = c("#1515FF", "#CECE0E", "#918D8D"),
              add = "reg.line", xlim = c(5, 14)) +
    stat_cor(aes(color = condition), method = "spearman", label.x = 8),
  nrow = 1
)

# Save the plot as SVG
C20H34O5_scaterplot 
C5H11NO2_scaterplot

ggsave(paste0("output/C20H34O5_scaterplot_season", year, ".svg"), 
       C20H34O5_scaterplot  , device = "svg")
ggsave(paste0("output/C20H34O5_scaterplot_season", year, ".png"), 
       C20H34O5_scaterplot  , device = "png")

ggsave(paste0("output/C5H11NO2_scaterplot_season", year, ".svg"), 
       C5H11NO2_scaterplot  , device = "svg")
ggsave(paste0("output/C5H11NO2_scaterplot_season", year, ".png"), 
       C5H11NO2_scaterplot  , device = "png")
