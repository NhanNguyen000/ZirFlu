# check sample with high cytokine value --------------------
k <- ZirFlu$log2_cytokine %>% filter(GRO....CXCL1.>0)
ZirFlu$log2_cytokine %>% filter(GRO....CXCL1.>0) %>% select(PatientID)
ZirFlu$log2_cytokine %>% filter(IFN..>0) %>% select(PatientID)
ZirFlu$log2_cytokine %>% filter(IFN...1 > 4) %>% select(PatientID)


ZirFlu$log2_abTiters_Feb2022

b <- cor(ZirFlu$log2_cytokine[4:37], ZirFlu$log2_abTiters_Feb2022$H1N1_d0)