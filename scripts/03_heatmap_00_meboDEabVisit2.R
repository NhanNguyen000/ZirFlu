rm(list = ls())

library(tidyverse)
library(ComplexHeatmap)

# functions used in the analysis -----------------------------------------------

# load the lm() model outcome --------------------------------------------------
# lm() model of antibody titer visit2 associated with metabolites
load("data/lmRes_meboAbVisit2.RData")

# DE metabolites correlated with antibody titer visit2
load("data/DEmebo_abVisit2.RData")

# heatmap for all metabolites -----------------------------------------------
lmStatistic <- lmRes %>% 
  lapply(function(x) x %>% select(targetVariable, statistic, Formula)) %>%
  imap(~cbind(.x, ab = .y))

lmStatistic <- lmRes %>% 
  imap(~cbind(.x, ab = .y)) %>% purrr::reduce(full_join) %>% 
  mutate(season = substring(ab, 1, 4)) %>%
  mutate(significance = ifelse(p.value < 0.05, TRUE, FALSE))

plotDat_wide <- lmStatistic %>%
  select(statistic, Formula, ab) %>% 
  pivot_wider(names_from = "ab", values_from = statistic) %>%
  column_to_rownames("Formula")

#plotDat_wide %>% as.matrix() %>% heatmap(Colv = NA, Rowv = NA)
#plotDat_wide %>% as.matrix() %>% Heatmap()
plotDat_wide %>% as.matrix() %>% Heatmap(cluster_columns = FALSE)

lmStatistic %>% 
  ggplot(aes(x = ab, y = Formula)) + 
  geom_tile(aes(fill = statistic)) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# heatmap for for selected DE metabolites --------------------------------------

## prepare input data for ploting heatmap  -------------------------------------------

# selected DE metabolites in season 2019 which had sig. associations with >= 2 out of 4 antibodies
DEmebo_2019 <- venn_meboAbVisit2_2019 %>% flatten() %>% unlist()
length(unique(DEmebo_2019))
selected_DEmebo_2019 <- DEmebo_2019[duplicated(DEmebo_2019)]

DEmebo_2020 <- venn_meboAbVisit2_2020 %>% flatten() %>% unlist()
length(unique(DEmebo_2020))

# only select metabolites show consistent trend (>=75% across strains and years)
consistented_metabolites <- lmStatistic %>% group_by(Formula) %>%
  summarise(postitive = sum(statistic > 0, na.rm = TRUE)) %>%
  full_join(lmStatistic %>% group_by(Formula) %>%
              summarise(negative = sum(statistic < 0, na.rm = TRUE))) %>%
  mutate(select = ifelse(postitive >= 8*0.75 | negative >= 8*0.75, TRUE, NA)) %>%
  filter(select == TRUE)

selected_lmStatistic <- lmStatistic %>% 
  filter(Formula %in% consistented_metabolites$Formula) %>%
  filter(Formula %in% selected_DEmebo_2019)

save(selected_lmStatistic, file = "data/DEmebo_abVisit2_consistAcrossStrainSeason.RData")

## plot heatmap per each season -------------------------------------------------------------
year <- "2019"
#year <- "2020"

plotDat <- selected_lmStatistic %>% filter(season == year) %>% 
  mutate(ab = gsub("2019_", "", ab),
         ab = gsub("_Visit2", "", ab)) %>%
  mutate(ab = factor(ab, levels = c("H1N1", "H3N2", "Bvictoria", "Byamagata"))) 

plotDat_wide <- plotDat %>% select(ab, statistic, Formula) %>%
  pivot_wider(names_from = "ab", values_from = statistic) %>% column_to_rownames("Formula")

heatmap_DEmebo <- plotDat %>%
  ggplot(aes(x = ab, y = Formula)) + 
  geom_tile(aes(fill = statistic)) +
  geom_text(aes(label = ifelse(significance == TRUE, "*", ""))) +
  scale_y_discrete(limits = plotDat$Formula[hclust(dist(plotDat_wide))$order])+
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        text = element_text(size = 12, face = "bold"))

# Save the plot as SVG
heatmap_DEmebo
ggsave(paste0("output/heatmap_DEmebo_season", year, ".svg"), 
       heatmap_DEmebo, device = "svg")
ggsave(paste0("output/heatmap_DEmebo_season", year, ".png"), 
       heatmap_DEmebo, device = "png")

## plot both season in 1 heatmap ------------------------------------------------------------
selected_lmStatistic %>%
  ggplot(aes(x = ab, y = Formula)) + 
  geom_tile(aes(fill = statistic)) +
  geom_text(aes(label = ifelse(significance == TRUE, "*", ""))) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
