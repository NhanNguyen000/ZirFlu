library(dendextend)

dend <- ZirFlu$abTiters_Feb2022 %>%
  dist() %>% hclust("ward.D2") %>% as.dendrogram()

dend <- ZirFlu$metabolite_dat %>%
  dist() %>% hclust("ward.D2") %>% as.dendrogram()
par(mar = c(3,30,1,20),  cex=0.7)
plot(dend, horiz = TRUE)
axis(side = 1, at = seq(0, 40000, 10000))
