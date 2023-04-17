library(caret)
library(tidyverse)

# function need for this code -------------------------------------------
get.cowplot_pca <- function(pca_dat, metadata) {
  cowplot::plot_grid(
    get.pca_plot(pca_dat, metadata, "condition"),
    get.pca_plot(pca_dat, metadata, "disease"),
    get.pca_plot(pca_dat, metadata, "category"),
    get.pca_plot(pca_dat, metadata, "time"),
    get.pca_plot(pca_dat, metadata, "sex"),
    get.pca_plot(pca_dat, metadata, "age_class")
  )
}

# make PCA for protein and metabolite data -------------------------------------------
metadata <- ZirFlu$donorSamples %>% full_join(ZirFlu$donorInfo) %>%
  mutate(disease = ifelse(condition == "healthy", "healthy", "cirrhosis")) %>%
  mutate(age_class = ifelse(age >60, ">60", "young"))

protein_pca <- predict(preProcess(x=ZirFlu$protein_dat, method = "knnImpute"),
                       ZirFlu$protein_dat) %>% prcomp()
summary(protein_pca)$importance[1:3, 1:5]
get.cowplot_pca(pca_dat = protein_pca, metadata)

metabolite_pca <- ZirFlu$metabolite_dat %>% prcomp()
summary(metabolite_pca)$importance[1:3, 1:5]
get.cowplot_pca(pca_dat = metabolite_pca, metadata)

# put select PCA plot together -------------------------------------
Fig1.C <- cowplot::plot_grid(
  get.pca_plot(protein_pca, metadata, "condition"),
  get.pca_plot(metabolite_pca, metadata, "condition"),
  nrow = 2, labels = c("protein", "metabolite"), 
  label_fontface = "plain", label_x = 0.5
)
par(mar = c(3,3,1,1))
plot(Fig1.C)