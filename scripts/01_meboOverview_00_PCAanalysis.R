rm(list = ls())

library(tidyverse)
library(caret)

# functions used in the analysis -----------------------------------------------

get.pca_plot <- function(pca, metadat, groupType) {
  library(tidyverse)
  library(caret)
  
  ggplot(data.frame(pca$x), aes(PC1, PC2, color = metadat[[groupType]])) +
    geom_point(alpha=.6, size = 2.5) +
    theme_classic() +
    theme(aspect.ratio = 1, plot.title = element_text(hjust = .5)) +
   # ggsci::scale_color_nejm() +
    labs(x = paste0("PC1 [", 
                    round(100*summary(pca)$importance["Proportion of Variance", "PC1"], 
                          digits = 2), 
                    "%]"),
         y = paste0("PC2 [", 
                    round(100*summary(pca)$importance["Proportion of Variance", "PC2"], 
                          digits = 2), 
                    "%]"), 
         color = groupType) +
    stat_ellipse()
}


# load the data list object------------------------------------------------------
load("processedData/ZirFlu.RData")

# make PCA plot for each year
year <- "2019" # or
year <- "2020"

# PCA for health condition -------------------------------------------------------
metadata_condition <- ZirFlu$donorSamples %>% 
  full_join(ZirFlu$donorInfo) %>% 
  drop_na() %>% filter(season == year)

metadata_condition %>% filter(probenID %in% rownames(ZirFlu$metabolite_dat)) %>% 
  count(season, condition, time)

mebo_pca_condition <- ZirFlu$metabolite_dat %>% 
  rownames_to_column("probenID") %>% 
  filter(probenID %in% metadata_condition$probenID) %>%
  column_to_rownames("probenID") %>%
  prcomp()
summary(mebo_pca_condition)$importance[1:3, 1:5]

# cowplot::plot_grid(
#   get.pca_plot(pca = mebo_pca_condition, metadat = metadata_condition, groupType = "sex"),
#   get.pca_plot(pca = mebo_pca_condition, metadat = metadata_condition, groupType = "time"),
#   nrow = 2
# )

pcaPlot_condition <- get.pca_plot(pca = mebo_pca_condition, 
             metadat = metadata_condition, 
             groupType = "condition") +
  scale_color_manual(values = c("Healthy" = "#0771C3", 
                                "Compensated cirrhosis" = "#ECC02C", 
                                "Decompensated cirrhosis" = "#A1A1A1")) +
  theme(legend.position = "top", text = element_text(size=14, face = "bold"))

# Save the plot as SVG
pcaPlot_condition 
ggsave(paste0("output/pcaPlot_condition_season", year, ".svg"), 
       pcaPlot_condition, device = "svg")
ggsave(paste0("output/pcaPlot_condition_season", year, ".png"), 
       pcaPlot_condition, device = "png")

# PCA for time visit -----------------------------------------------------------------------------
metadata_visit <- ZirFlu$donorSamples %>% 
  full_join(ZirFlu$donorInfo) %>% 
  drop_na() %>% 
  filter(season == year) %>%
  filter(time != "Visit1") %>%
  mutate(time = ifelse(time == "Visit2", "Visit 2", time))

metadata_visit %>% filter(probenID %in% rownames(ZirFlu$metabolite_dat)) %>% 
  count(season, condition, time)

mebo_pca_visit <- ZirFlu$metabolite_dat %>% 
  rownames_to_column("probenID") %>% 
  filter(probenID %in% metadata_visit$probenID) %>%
  column_to_rownames("probenID") %>%
  prcomp()
summary(mebo_pca_visit)$importance[1:3, 1:5]

pcaPlot_visit <- get.pca_plot(pca = mebo_pca_visit, 
             metadat = metadata_visit, 
             groupType = "time") +
  theme(legend.position = "top", text = element_text(size=14, face = "bold"))

# Save the plot as SVG
pcaPlot_visit
ggsave(paste0("output/pcaPlot_visit_season", year, ".svg"), 
       pcaPlot_visit, device = "svg")
ggsave(paste0("output/pcaPlot_visit_season", year, ".png"), 
       pcaPlot_visit, device = "png")
