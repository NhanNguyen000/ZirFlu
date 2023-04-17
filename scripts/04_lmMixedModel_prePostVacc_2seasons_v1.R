library(tidyverse)
library(limma)

## protein data ------------------------
inputDat_protein <- ZirFlu$donorSamples %>% 
  full_join(ZirFlu$donorInfo) %>% right_join(ZirFlu$HAItiter) %>% 
  left_join(ZirFlu$protein_dat %>% rownames_to_column("probenID"))

inputDat_protein_NAcount <- inputDat_protein %>% group_by(time, disease) %>%
  summarise(across(everything(), ~ sum(is.na(.))))
inputDat_protein %>% group_by(season, time, disease) %>% count()

k <- inputDat_protein_NAcount[,c(1, 2)] %>%
  cbind(inputDat_protein_NAcount[,-c(1, 2)]/inputDat_sampleCount$n)

k1 <- k %>% filter(time == "T2") %>% 
  summarise(across(starts_with("OID"), ~ .x >=0.3))
a <- names(which(colSums(k1[, -c(1)]) >0))

k1 <- k %>% filter(time != "T3") %>%
  summarise(across(starts_with("OID"), ~ .x >=0.3)) %>% t()

k2 <- k %>% filter(time != "T2")

# per season 
input_protein <- list() # use missing data
input_protein_complete <- list() # use complete sample have both pre and post-vaccine measurement.
for(year in c("2019", "2020")) {
  input_protein[[year]]$T1vsT2 <- inputDat_protein %>% 
    filter(season == year) %>% filter(time %in% c("T1", "T2"))
  input_protein[[year]]$T1vsT3 <- inputDat_protein %>% 
    filter(season == year) %>% filter(time %in% c("T1", "T3"))
  
  input_protein_complete[[year]]$T1vsT2 <- inputDat_protein %>% 
    filter(season == year) %>% filter(time %in% c("T1", "T2")) %>%
    add_count(patientID) %>% filter(n == 2) %>% select(-n)
  input_protein_complete[[year]]$T1vsT3 <- inputDat_protein %>% 
    filter(season == year) %>% filter(time %in% c("T1", "T3")) %>%
    add_count(patientID) %>% filter(n == 2) %>% select(-n)
}

# check NAs in protein data
a <- input_protein$`2019`$T1vsT2
tidyv
# continute
lmRes_prePost_protein <- list()
for (season in names(input_protein)) {
  for (times in names(input_protein[[season]])) {
    metaDat <- input_protein[[season]][[times]] %>%
      select(patientID, probenID, time, sex, age, disease, category)
    
    designTable <- model.matrix(~ time + sex + age + disease, metaDat)
    
    proteinDat <- input_protein[[season]][[times]] %>% 
      select(-c(1:26)) %>% t()
    
    dupcor <- duplicateCorrelation(object = proteinDat, design = designTable,
                                   block = metaDat$patientID)
    
    lmRes_prePost_protein[[season]][[times]] <- lmFit(proteinDat, 
                                                      design = designTable, 
                                                      block = metaDat$patientID,
                                                      correlation = dupcor$consensus.correlation) %>% eBayes() 
  }
}

lmRes_prePost_protein$`2019`$T1vsT2 %>% decideTests() %>% summary()
lmRes_prePost_protein$`2019`$T1vsT3 %>% decideTests() %>% summary()
lmRes_prePost_protein$`2020`$T1vsT2 %>% decideTests() %>% summary()
lmRes_prePost_protein$`2020`$T1vsT3 %>% decideTests() %>% summary()

resPrePost_protein <- list()
for(season in names(lmRes_prePost_protein)) {
  for(times in names(lmRes_prePost_protein[[season]])) {
    resTable <- lmRes_prePost_protein[[season]][[times]] %>% 
      topTable(coef = 2, sort.by = 'none', adjust.method = 'BH', number = Inf)
    resPrePost_protein[[season]][[times]]$resTable <- resTable
    resPrePost_protein[[season]][[times]]$pval <- resTable %>% filter(P.Value < 0.05)
    resPrePost_protein[[season]][[times]]$padj <- resTable %>% filter(adj.P.Val < 0.05)
  }
}

a<-intersect(rownames(resPrePost_protein$`2019`$T1vsT2$padj),
          rownames(resPrePost_protein$`2019`$T1vsT3$padj)) # 21 proteins

A <- get.proteinAnnot(ZirFlu$protein_annot,
                      rownames(resPrePost_protein$`2019`$T1vsT3$padj)) %>%
  full_join(resPrePost_protein$`2019`$T1vsT3$padj %>% rownames_to_column("OlinkID"))

intersect(a, rownames(resPrePost_protein$`2020`$T1vsT2$pval)) # no overlap
intersect(a, rownames(resPrePost_protein$`2020`$T1vsT3$pval)) # 2 overlap

a<- get.proteinAnnot(ZirFlu$protein_annot, rownames(resPrePost_protein$`2019`$T1vsT2$padj)) 
b<- get.proteinAnnot(ZirFlu$protein_annot, rownames(resPrePost_protein$`2019`$T1vsT3$padj)) 
selected_proteins <- unique(c(rownames(resPrePost_protein$`2019`$T1vsT2$padj),
                            rownames(resPrePost_protein$`2019`$T1vsT3$padj)))

rownames(resPrePost_protein$`2019`$T1vsT3$pval)
rownames(resPrePost_protein$`2019`$T1vsT3$padj)
# need to check later ---------------------------
get.proteinAnnot(ZirFlu$protein_annot, rownames(outcome$T3$padj))
get.proteinAnnot(ZirFlu$protein_annot, rownames(outcome$T4$padj))

# other code
topTable(lmRes, coef = 2)
topTable(lmRes, coef = "timeT4")

topTable(lmRes, coef = 5)
topTable(lmRes, coef = "diseasehealthy")

volcanoplot(lmRes, coef = "timeT4", highlight = 10)
abline(h = -log10(0.05))

plotMD(lmRes, column = 2)
qqt(lmRes$t[,2], df = lmRes$df.residual + lmRes$df.prior)
abline(0,1)

## metabolite data ------------------------
inputDat_metabolite <- ZirFlu$donorSamples %>% 
  full_join(ZirFlu$donorInfo) %>% right_join(ZirFlu$HAItiter) %>% 
  left_join(ZirFlu$metabolite_dat %>% rownames_to_column("probenID"))

inputDat_metabolite %>% group_by(season, time, disease) %>% count()

# per season 
input_metabolite <- list()
for(year in c("2019", "2020")) {
  input_metabolite[[year]]$T1vsT2 <- inputDat_metabolite%>% 
    filter(season == year) %>% filter(time %in% c("T1", "T2"))
  input_metabolite[[year]]$T1vsT3 <- inputDat_metabolite %>% 
    filter(season == year) %>% filter(time %in% c("T1", "T3"))
}

# continute
lmRes_prePost_metabolite <- list()
for (season in names(input_metabolite)) {
  for (times in names(input_metabolite[[season]])) {
    metaDat <- input_metabolite[[season]][[times]] %>%
      select(patientID, probenID, time, sex, age, disease, category)
    
    designTable <- model.matrix(~ time + sex + age + disease, metaDat)
    
    metaboliteDat <- input_metabolite[[season]][[times]] %>% 
      select(-c(1:26)) %>% t()
    
    dupcor <- duplicateCorrelation(object = metaboliteDat, design = designTable,
                                   block = metaDat$patientID)
    
    lmRes_prePost_metabolite[[season]][[times]] <- lmFit(metaboliteDat, 
                                                      design = designTable, 
                                                      block = metaDat$patientID,
                                                      correlation = dupcor$consensus.correlation) %>% eBayes() 
  }
}

lmRes_prePost_metabolite$`2019`$T1vsT2 %>% decideTests() %>% summary()
lmRes_prePost_metabolite$`2019`$T1vsT3 %>% decideTests() %>% summary()
lmRes_prePost_metabolite$`2020`$T1vsT2 %>% decideTests() %>% summary()
lmRes_prePost_metabolite$`2020`$T1vsT3 %>% decideTests() %>% summary()

resPrePost_metabolite <- list()
for(season in names(lmRes_prePost_metabolite)) {
  for(times in names(lmRes_prePost_metabolite[[season]])) {
    resTable <- lmRes_prePost_metabolite[[season]][[times]] %>% 
      topTable(coef = 2, sort.by = 'none', adjust.method = 'BH', number = Inf)
    resPrePost_metabolite[[season]][[times]]$resTable <- resTable
    resPrePost_metabolite[[season]][[times]]$pval <- resTable %>% filter(P.Value < 0.05)
    resPrePost_metabolite[[season]][[times]]$padj <- resTable %>% filter(adj.P.Val < 0.05)
  }
}

