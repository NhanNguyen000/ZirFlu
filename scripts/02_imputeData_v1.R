## protein - impute missing value  ---------------------------------------------
library(naniar)
library(visdat)

# check the mising in the sample
check.NAs(ZirFlu$protein_dat)
# prop_complete(ZirFlu$protein_dat)
# prop_miss(ZirFlu$protein_dat)
# miss_var_table(ZirFlu$protein_dat)
# miss_case_table(ZirFlu$protein_dat)

# remove proteins have more than > 10% NAs
ZirFlu$protein_dat2 <- ZirFlu$protein_dat %>% 
  select(which(colSums(is.na(ZirFlu$protein_dat))/nrow(ZirFlu$protein_dat) < 0.1))
check.NAs(ZirFlu$protein_dat2)

## Impute missing value with median (MICE package only supports mean & regression) -----------------------------
library(mice)
# Impute the missing value without considering longitudinal aspect (wide data format)
ZirFlu$protein_imputedDat <- ZirFlu$protein_dat %>% 
  mutate_if(is.numeric, function(x) ifelse(is.na(x), median(x, na.rm = TRUE), x))

# # Impute the missing value without considering longitudinal aspect
# option "pmm" is the default - apply for all data type - take a lot of time
# option "norm" (Bayesian linear regression, recommended method) - could not run because have a number of unbalanced factors
# option "norm.nob" (linear regression, ignoring model error) - able to runs
# option "norm.boot" (linear regression, using bootstrap)
# option "norm.predict" (linear regression, predicted values)
# Impute the missing value without considering longitudinal aspect - with different methods

init <- mice(ZirFlu$protein_dat, maxit = 0)
predM <- init$predictorMatrix
imputed_dat <- list()
for (selected_method in c("pmm", "norm.nob", "norm.boot", "norm.predict")) {
  meth <- rep(selected_method, length(init$method))
  imputed <- mice(ZirFlu$protein_dat, method = meth, predictorMatrix = predM, m = 5)
  imputed_dat[[selected_method]] <- complete(imputed)
}

init <- mice(ZirFlu$protein_dat, maxit = 0)
predM <- init$predictorMatrix
imputed_dat <- list()
for (selected_method in c("norm.predict")) {
  meth <- rep(selected_method, length(init$method))
  imputed <- mice(ZirFlu$protein_dat, method = meth, predictorMatrix = predM, m = 5)
  imputed_dat[[selected_method]] <- complete(imputed)
}
save(imputed_dat, "20220909_ZirFlu_imputedProteinsDat_NN.RData")

## Protein - long format (same donor with measurement at T1, T2, T3 - longitudinal data) ----------------------------
### not remove sample have more than >50% NAs ----------------------------
proteinImputeDat_time <- list()
for (time_point in c("T1", "T3", "T4")) {
  proteinImputeDat_time[[time_point]] <- ZirFlu$protein_dat %>% 
    rownames_to_column("probenID") %>% 
    filter(probenID %in% ZirFlu$metadata2$probenID[which(ZirFlu$metadata2$Time == time_point)]) %>%
    mutate_if(is.numeric, function(x) ifelse(is.na(x), median(x, na.rm = TRUE), x)) %>%
    column_to_rownames("probenID")
}

dat_temp <- proteinImputeDat_time %>% reduce(rbind)
ZirFlu$proteinImputeDat <- dat_temp[rownames(ZirFlu$protein_dat),]
identical(rownames(ZirFlu$proteinImputeDat), rownames(ZirFlu$protein_dat))

ZirFlu$proteinImputeDat_time <- proteinImputeDat_time

### remove sample have more than >50% NAs ----------------------------
proteinImputeDat_time2 <- list()
for (time_point in c("T1", "T3", "T4")) {
  dat_temp <- ZirFlu$protein_dat %>% 
    rownames_to_column("probenID") %>% 
    filter(probenID %in% ZirFlu$metadata2$probenID[which(ZirFlu$metadata2$Time == time_point)]) %>%
    column_to_rownames("probenID")
  
  proteinImputeDat_time2[[time_point]] <- dat_temp %>%
    select(which(colSums(is.na(dat_temp))/nrow(dat_temp) < 0.1)) %>%
    mutate_if(is.numeric, function(x) ifelse(is.na(x), median(x, na.rm = TRUE), x))
}

ZirFlu$proteinImputeDat_time2 <- proteinImputeDat_time2
