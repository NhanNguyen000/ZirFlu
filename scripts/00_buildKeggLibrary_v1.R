library("KEGGREST")
library(tidyverse)

cpdpathway <- keggLink("cpd", "pathway") %>% enframe %>%
  rename("mapId" = "name", "cpdId" = "value")

listpathway <- keggList("pathway") %>% enframe %>%
  rename("mapId" = "name", "pathway" = "value")

# cpdpathway_combine <- cpdpathway %>% full_join(listpathway) %>%
#   group_by(cpdId) %>% summarise(across(everything(), str_c, collapse=",")) 

cpdpathway_combine <- cpdpathway %>% full_join(listpathway)

keggID_library <-ZirFlu$metabolite_annot %>% 
  slice(-grep("HMDB|CHEBI", CompoundID)) %>%
  inner_join(cpdpathway_combine %>% mutate(cpdId = substring(cpdId, 5, 10)), 
            by = c("CompoundID" = "cpdId")) %>%
  select(pathway, CompoundID) %>% group_by(pathway) %>%
  summarise(across(everything(), str_c, collapse = "; "))
# 293 pathway in total, 375 metabolite Idx
write.csv(keggID_library, file = "temp/20221013_keggID_library.csv", 
          row.names = FALSE)

keggID_library2 <-ZirFlu$metabolite_annot %>% 
  slice(-grep("HMDB|CHEBI", CompoundID)) %>%
  inner_join(cpdpathway_combine %>% mutate(cpdId = substring(cpdId, 5, 10)), 
             by = c("CompoundID" = "cpdId")) %>%
  select(pathway, CompoundName) %>% group_by(pathway) %>%
  summarise(across(everything(), str_c, collapse = "; "))
write.csv(keggID_library2, file = "temp/20221013_keggID_library2.csv", 
          row.names = FALSE)

# test with Jianbo code
fgid <- read.xlsx2("temp/300BCG.all.HMDB.KEGG1.xlsx", 1, header=TRUE) 
TEST <- fgid %>% # 386, but I only have 337
  mutate(KEGG= paste0("cpd:", KEGG)) %>%
  inner_join(cpdpathway_combine, 
             by = c("KEGG" = "cpdId")) %>% 
  select(pathway, Match) %>% group_by(pathway) %>%
  summarise(across(everything(), str_c, collapse = "; "))
# note, both of his library are the same, just 1 library double the name

# Jianbo code: 
#BiocManager::install("KEGGREST",force = TRUE)
library(KEGGREST)
cpdpathway = keggLink("cpd", "pathway")
class(cpdpathway)
head(cpdpathway,n=5)
length(cpdpathway)

listpathway = keggList("pathway")
head(listpathway,n=5)
listpathway1 = as.data.frame(listpathway)
head(listpathway1,n=5)
row.names(listpathway1)
listpathway1$mapid = row.names(listpathway1)
row.names(listpathway1) <- NULL
#save(listpathway1,file = "mapid.pathwayname.rda")
#row.names(listpathway1)
#先把mapid转换成pathway names
#再按cpd id合并，即一个cpd id对应多个 pathway name
bb <- as.matrix(cpdpathway)
head(bb)
b1 <- row.names(bb)
bb1 <- cbind(b1,bb)
bb2 = as.data.frame(bb1)
head(bb2,n=3)
row.names(bb2) <- NULL
names(bb2)[1] <- "mapid"
names(bb2)[2] <- "cpdid"
dim(bb2)
#save(bb2,file = "mapid.cpdid.rda")

merge1 <- merge(bb2, listpathway1, by = "mapid", alL = FALSE)
dim(merge1)
head(merge1,n=3)

#合并相同行名
library(dplyr)
library(stringr)
bbb <- merge1 %>% 
  group_by(cpdid) %>% 
  summarise(across(everything(), str_c, collapse=",")) 
head(bbb,n=3)
#save(bbb,file = "all.cpdid.mapid.pathway.rda")
library("xlsx")
fgid <- read.xlsx2("300BCG.all.HMDB.KEGG1.xlsx", 1, header=TRUE)
head(fgid)
library(dplyr)
fgid = fgid %>% distinct(Match,KEGG, .keep_all = TRUE)
dim(fgid)
fgid$cpd = c(rep("cpd",1760))
fgid$cpdid <- paste(fgid$cpd,fgid$KEGG, sep = ":")
xx  = c("cpd","KEGG")
fgid1 = fgid[,!names(fgid) %in% xx]
head(fgid1)
library(dplyr)
#fgid2 <- fgid1 %>% distinct(cpdid, .keep_all = T)
dim(bbb)
dim(fgid1)
merge111 <- merge(fgid1, bbb, by = "cpdid", alL = FALSE)
dim(merge111)
head(merge111,n=3)
merge111$cpdid
fgid1$cpdid
ad1 <- fgid1$cpdid
bd1 <- bbb$cpdid
bba <- Reduce(intersect, list(ad1,bd1))
length(bba)
aas <- setdiff(ad1, bd1)

xx  = c("mapid")
merge112 = merge111[,!names(merge111) %in% xx]
#merge112$listpathway
#拆分一行变多行
library(tidyr)
dff2 <-merge112 %>% as_tibble() %>% 
  separate_rows(listpathway, sep = ",")
head(dff2,n=5)

dff3 <- dff2 %>% distinct(listpathway, .keep_all = T)

cccc <- dff2 %>% 
  group_by(listpathway) %>% 
  summarise(across(everything(), str_c, collapse="; ")) 
head(cccc,n=3)
#cccc$Match
dim(cccc)
xx  = c("cpdid")
cccc1 = cccc[,!names(cccc) %in% xx]
write.csv(cccc1,file = "300BCG.kegglib.csv",row.names = FALSE)

ddff <- read.csv("300BCG.kegglib.csv",header=TRUE,sep=",",stringsAsFactors = FALSE)
head(ddff,n=3)
dim(ddff)

ddff1 <-ddff %>% as_tibble() %>% 
  separate_rows(Match, sep = "; ")
head(ddff1,n=3)
dim(ddff1)
ddff2 <- unique(ddff1)
dim(ddff2)

ddff3 <- ddff2 %>% 
  group_by(listpathway) %>% 
  summarise(across(everything(), str_c, collapse="; ")) 

dim(ddff3)
write.csv(ddff3,file = "300BCG.kegglib1.csv",row.names = FALSE)