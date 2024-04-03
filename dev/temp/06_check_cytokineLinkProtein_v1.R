cytokineName <- "CXCL1"
cytokineName <- c("CXCL1", )
annot_temp <- ZirFlu$protein_annot %>% filter(Assay %in% cytokineName)

proteinCytokine_dat <- ZirFlu$protein_dat %>% 
  select(annot_temp$OlinkID) %>%
  rownames_to_column("probenID") %>% 
  mutate(probenID = as.numeric(probenID)) %>%
  left_join(ZirFlu$metadata %>% select(Patient.ID, probenID)) %>%
  left_join(ZirFlu$log2_cytokine, by = c("Patient.ID" = "PatientID"))

cor(proteinCytokine_dat$OID20762, proteinCytokine_dat$GRO....CXCL1.,
    use = 'pairwise.complete.obs')

proteinCytokine_dat %>% 
  select(probenID, OID20762,GRO....CXCL1.) %>% distinct() %>%
  ggplot(aes(OID20762,GRO....CXCL1.)) + geom_point() +
  xlab("OID20762 (CXCL1)") + ylab("Log2 of GRO....CXCL1.")

selected_proteinDat_abTiters <- proteinCytokine_dat %>% 
  left_join(ZirFlu$log2_abTiters_Feb2022, 
            by = c("Patient.ID" = "PatientID", "Condition"))

cowplot::plot_grid(
  get.scatterPlot2(selected_proteinDat_abTiters %>% filter(Measure == "d0"),
                   x_name = 'H1N1_d0', y_name = "GRO....CXCL1.",
                   color_group = "Condition"),
  get.scatterPlot2(selected_proteinDat_abTiters %>% filter(Measure == "d21-35"),
                   x_name = 'H1N1_d21.35', y_name = "GRO....CXCL1.",
                   color_group = "Condition"),
  get.scatterPlot2(selected_proteinDat_abTiters %>% filter(Measure == "d0"),
                   x_name = 'H3N2_d0', y_name = "GRO....CXCL1.",
                   color_group = "Condition"),
  get.scatterPlot2(selected_proteinDat_abTiters %>% filter(Measure == "d21-35"),
                   x_name = 'H3N2_d21.35', y_name = "GRO....CXCL1.",
                   color_group = "Condition"),
  get.scatterPlot2(selected_proteinDat_abTiters %>% filter(Measure == "d0"), 
                   x_name = 'Bvictoria.Maryland_d0', y_name = "GRO....CXCL1.",
                   color_group = "Condition"),
  get.scatterPlot2(selected_proteinDat_abTiters %>% filter(Measure == "d21-35"), 
                   x_name = 'Bvictoria.Maryland_d21.35', y_name = "GRO....CXCL1.",
                   color_group = "Condition"),
  get.scatterPlot2(selected_proteinDat_abTiters %>% filter(Measure == "d0"),
                   x_name = 'Byamagata.Phuket_d0', y_name = "GRO....CXCL1.",
                   color_group = "Condition"),
  get.scatterPlot2(selected_proteinDat_abTiters %>% filter(Measure == "d21-35"),
                   x_name = 'Byamagata.Phuket_d21.35', y_name = "GRO....CXCL1.",
                   color_group = "Condition"),
  nrow=2
)


