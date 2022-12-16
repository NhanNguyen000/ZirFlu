library(tidyverse)
library(limma)

get.metaDat <- function(cohort, times) {
  dat <- cohort$donorInfo %>% full_join(cohort$donorSamples) %>%
    select(patientID, probenID, time, sex, age, condition, category) %>%
    mutate(disease = ifelse(condition == "healthy", "healthy", "cirrhosis")) %>%
    slice(which(time %in% times))
  return(dat)
}

get.proteinVal <- function(cohort, selectSamples) {
  proteinVal <- cohort$protein_dat %>% rownames_to_column("probenID") %>%
    mutate(probenID = as.numeric(probenID)) %>%
    slice(which(probenID %in% selectSamples)) %>% 
    column_to_rownames("probenID") %>% t()
  return(proteinVal)
}

metaDat <- get.metaDat(ZirFlu, c("T1", "T4")) # use missing data
metaDat2 <- get.metaDat(ZirFlu, c("T1", "T4")) %>% 
  add_count(patientID) %>% filter(n == 2) %>% select(-n) # use complete sample have both pre and post-vaccine measurement.

proteinDat <- get.proteinVal(ZirFlu, metaDat$probenID)

designVal <- metaDat %>% arrange(match(probenID, colnames(proteinDat)))
identical(as.numeric(colnames(proteinDat)), as.numeric(designVal$probenID))

designTable <- model.matrix(~ time + sex + age + disease, designVal)
dupcor <- duplicateCorrelation(object = proteinDat, design = designTable,
                               block = designVal$patientID)
lmRes <- lmFit(proteinDat, design = designTable,
               block = designVal$patientID,
               correlation = dupcor$consensus.correlation) %>% eBayes() 

lmRes %>% decideTests() %>% summary()
resTable <- lmRes %>% 
  topTable(coef = paste0("time", postTime), sort.by = 'none', adjust.method = 'BH', number = Inf)

padj <- resTable %>% filter(adj.P.Val < 0.05)
pval <- resTable %>% filter(P.Value < 0.05)

# could use another oder?: protein = time + patient + disease + sex + age
# complete cases - sould consider Patient as ramdom effect?

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


