## antibody (ab) titers --------------------------------
log2_ab <- ZirFlu$abTiters_Feb2022 %>% 
  mutate(across(where(is.numeric), ~log2(.x))) %>%
  mutate(across(where(is.numeric), ~replace_na(.,0))) %>%
  mutate(across(where(is.numeric), ~ifelse(is.infinite(.x), 0, .x)))

cowplot::plot_grid(
  get.comparMeasures_plot(log2_dat = log2_ab, 
                          x_name = "H1N1_d0", y_name = "H1N1_d21.35", 
                          color_group = "Condition", abline_val = 41),
  get.comparMeasures_plot(log2_dat = log2_ab, 
                          x_name = "H3N2_d0", y_name = "H3N2_d21.35", 
                          color_group = "Condition", abline_val = 10),
  get.comparMeasures_plot(log2_dat = log2_ab, 
                          x_name = "Bvictoria.Maryland_d0", y_name = "Bvictoria.Maryland_d21.35", 
                          color_group = "Condition", abline_val = 10),
  get.comparMeasures_plot(log2_dat = log2_ab, 
                          x_name = "Byamagata.Phuket_d0", y_name = "Byamagata.Phuket_d21.35", 
                          color_group = "Condition", abline_val = 10)
)
cowplot::plot_grid(
  get.comparMeasures_plot(log2_dat = log2_ab, 
                          x_name = "H1N1_d0", y_name = "H1N1_d21.35", 
                          color_group = "Condition", abline_val = 40),
  get.comparMeasures_plot(log2_dat = log2_ab, 
                          x_name = "H3N2_d0", y_name = "H3N2_d21.35", 
                          color_group = "Condition", abline_val = 40),
  get.comparMeasures_plot(log2_dat = log2_ab, 
                          x_name = "Bvictoria.Maryland_d0", y_name = "Bvictoria.Maryland_d21.35", 
                          color_group = "Condition", abline_val = 40),
  get.comparMeasures_plot(log2_dat = log2_ab, 
                          x_name = "Byamagata.Phuket_d0", y_name = "Byamagata.Phuket_d21.35", 
                          color_group = "Condition", abline_val = 40)
)

log2_ab2 <- log2_ab %>% full_join(ZirFlu$metadata, by = c("PatientID" = "Patient.ID"))

cowplot::plot_grid(
  get.scatter_plot(log2_ab2, 
                   x_name = "H1N1_d0", y_name = "H1N1_d21.35",
                   color_group = "condition", shape_group = "sex"),
  get.scatter_plot(log2_ab2, 
                   x_name = "H1N1_d0", y_name = "H1N1_d21.35",
                   color_group = "condition", shape_group = "sex"),
  get.scatter_plot(log2_ab2, 
                   x_name = "Bvictoria.Maryland_d0", y_name = "Bvictoria.Maryland_d21.35", 
                   color_group = "condition", shape_group = "sex"),
  get.scatter_plot(log2_ab2, 
                   x_name = "Byamagata.Phuket_d0", y_name = "Byamagata.Phuket_d21.35", 
                   color_group = "condition", shape_group = "sex")
)


ZirFlu$log2_abTiters_Feb2022 <- ZirFlu$abTiters_Feb2022 %>% 
  get.log2() %>%
  mutate(Condition = ifelse(Condition == "Healthy", "Healthy", "Cirrhosis")) %>%
  mutate(ab_H1N1 = H1N1_d21.35- H1N1_d0,
         ab_H3N2 = H3N2_d21.35- H3N2_d0,
         ab_Bvictoria.Maryland = Bvictoria.Maryland_d21.35 - Bvictoria.Maryland_d0,
         ab_Byamagata.Phuket = Byamagata.Phuket_d21.35 - Byamagata.Phuket_d0)

cowplot::plot_grid(
  get.comparMeasures_plot(log2_dat = ZirFlu$log2_abTiters_Feb2022, 
                          x_name = "H1N1_d0", y_name = "H1N1_d21.35", 
                          color_group = "Condition", abline_val = 10),
  get.comparMeasures_plot(log2_dat = ZirFlu$log2_abTiters_Feb2022, 
                          x_name = "H3N2_d0", y_name = "H3N2_d21.35", 
                          color_group = "Condition", abline_val = 10),
  get.comparMeasures_plot(log2_dat = ZirFlu$log2_abTiters_Feb2022, 
                          x_name = "Bvictoria.Maryland_d0", y_name = "Bvictoria.Maryland_d21.35", 
                          color_group = "Condition", abline_val = 10),
  get.comparMeasures_plot(log2_dat = ZirFlu$log2_abTiters_Feb2022, 
                          x_name = "Byamagata.Phuket_d0", y_name = "Byamagata.Phuket_d21.35", 
                          color_group = "Condition", abline_val = 10)
)