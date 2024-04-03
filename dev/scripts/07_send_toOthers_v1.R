# send protein overvew to Valerie (2023.01.06) ----------------------
library("tidyverse")
library("xlsx")

protein_overview <- ZirFlu$protein_annot %>% 
  select(Assay, UniProt) %>% rename("Protein_name" = Assay)


write.xlsx(as.data.frame(protein_overview), file = "sendOut/20230106_proteinOverview_forValerie_NhanNguyen.xlsx",
           sheetName = "protein_overview", row.names = FALSE, append = FALSE)

# send beta-analine (C3H7NO2) and tyrosine (C9H11NO3) to Saumya ------------------
View(ZirFlu$metabolite_annot)

meta.Ids <- ZirFlu$metabolite_annot %>% 
  filter(Formula %in% c("C3H7NO2", "C9H11NO3")) %>%
  select(ionIdx, Formula) %>% distinct()

meta.Dat <- ZirFlu$metabolite_dat %>% 
  select(as.character(meta.Ids$ionIdx)) %>% 
  rename("C3H7NO2" = "25", "C9H11NO3" = "300") %>%
  rownames_to_column("probenID")

metadat <- ZirFlu$donorInfo %>% 
  full_join(ZirFlu$donorSamples) %>%
  full_join(ZirFlu$HAItiter) %>%
  full_join(meta.Dat) %>% drop_na(category) %>%
  mutate(category = factor(category, levels = c("NR", "Other", "TR")))

library(ggpubr)
metabolite <- "C3H7NO2"
metabolite <- "C9H11NO3"
ggboxplot(data = metadat, 
          x = "category", y = metabolite, 
          palette = "jco", add = "jitter") + 
  facet_wrap(~disease + season)

my_comparisons <- list(c("Other", "NR"), c("TR", "NR"), c("Other", "TR"))
ggboxplot(data = metadat, 
          x = "category", y = metabolite, 
          palette = "jco", add = "jitter") + 
  facet_wrap(~disease) + 
  stat_compare_means(comparisons = my_comparisons, method = "t.test")

