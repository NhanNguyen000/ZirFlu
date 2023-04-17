IL15protein_abTiters_Feb <- ZirFlu$protein_dat %>% select(OID20562) %>% 
  rownames_to_column("probenID") %>% 
  mutate(probenID = as.numeric(probenID)) %>%
  left_join(ZirFlu$metadata2 %>% 
              select(Patient.ID, probenID, Time, condition, category)) %>%
  filter(Time == "T1") %>%
  left_join(ZirFlu$log2_abTiters_Feb2022, by = c("Patient.ID" = "PatientID"))

IL15protein_abTiters_Aug <- ZirFlu$protein_dat %>% select(OID20562) %>% 
  rownames_to_column("probenID") %>% 
  mutate(probenID = as.numeric(probenID)) %>%
  left_join(ZirFlu$metadata2 %>% 
              select(Patient.ID, probenID, Time, condition, category)) %>%
  filter(Time == "T1") %>%
  left_join(ZirFlu$log2_abTiters_Aug2022, by = c("Patient.ID" = "PatientID")) %>%
  rename("IL-15 (OID20562)" = OID20562)

cowplot::plot_grid(
  get.scatterPlot3(IL15protein_abTiters_Feb, x_name = "H1N1_titerFC", 
                   y_name = "OID20562", color_group = "Category"),
  get.scatterPlot3(IL15protein_abTiters_Feb, x_name = "H3N2_titerFC", 
                   y_name = "OID20562", color_group = "Category"),
  get.scatterPlot3(IL15protein_abTiters_Feb, x_name = "Bvictoria.Maryland_titerFC", 
                   y_name = "OID20562", color_group = "Category"),
  get.scatterPlot3(IL15protein_abTiters_Feb, x_name = "Byamagata.Phuket_titerFC", 
                   y_name = "OID20562", color_group = "Category"),
  get.scatterPlot3(IL15protein_abTiters_Aug, x_name = "H1N1_titerFC", 
                   y_name = "OID20562", color_group = "Category"),
  get.scatterPlot3(IL15protein_abTiters_Aug, x_name = "H3N2_titerFC", 
                   y_name = "OID20562", color_group = "Category"),
  get.scatterPlot3(IL15protein_abTiters_Aug, x_name = "Bvictoria.Maryland_titerFC", 
                   y_name = "OID20562", color_group = "Category"),
  get.scatterPlot3(IL15protein_abTiters_Aug, x_name = "Byamagata.Phuket_titerFC", 
                   y_name = "OID20562", color_group = "Category"),
  nrow = 2
)

cowplot::plot_grid(
  get.scatterPlot4(IL15protein_abTiters_Aug, 
                   x_name = "H1N1_titerFC", 
                   y_name = "IL-15 (OID20562)", color_group = "condition"),
  get.scatterPlot4(IL15protein_abTiters_Aug, 
                   x_name = "H3N2_titerFC", 
                   y_name = "IL-15 (OID20562)", color_group = "condition"),
  get.scatterPlot4(IL15protein_abTiters_Aug, 
                   x_name = "Bvictoria.Maryland_titerFC", 
                   y_name = "IL-15 (OID20562)", color_group = "condition"),
  get.scatterPlot3(IL15protein_abTiters_Aug %>% filter(Condition == "Cirrhosis"), 
                   x_name = "H1N1_titerFC", 
                   y_name = "IL-15 (OID20562)", color_group = "Category"),
  get.scatterPlot3(IL15protein_abTiters_Aug %>% filter(Condition == "Cirrhosis"), 
                   x_name = "H3N2_titerFC", 
                   y_name = "IL-15 (OID20562)", color_group = "Category"),
  get.scatterPlot3(IL15protein_abTiters_Aug %>% filter(Condition == "Cirrhosis"), 
                   x_name = "Bvictoria.Maryland_titerFC", 
                   y_name = "IL-15 (OID20562)", color_group = "Category"),
  get.scatterPlot3(IL15protein_abTiters_Aug %>% filter(Condition == "Healthy"), 
                   x_name = "H1N1_titerFC", 
                   y_name = "IL-15 (OID20562)", color_group = "Category"),
  get.scatterPlot3(IL15protein_abTiters_Aug %>% filter(Condition == "Healthy"), 
                   x_name = "H3N2_titerFC", 
                   y_name = "IL-15 (OID20562)", color_group = "Category"),
  get.scatterPlot3(IL15protein_abTiters_Aug %>% filter(Condition == "Healthy"), 
                   x_name = "Bvictoria.Maryland_titerFC", 
                   y_name = "IL-15 (OID20562)", color_group = "Category"),
  nrow = 3
)
