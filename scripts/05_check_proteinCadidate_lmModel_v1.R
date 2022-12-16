inputDat <- ZirFlu$donorInfo %>% full_join(ZirFlu$donorSamples) %>%
  select(patientID, probenID, sex, age, condition, category, time) %>%
  full_join(ZirFlu$HAItiter_Feb2022 %>% select(-condition)) %>%
  right_join(ZirFlu$protein_dat %>% rownames_to_column("probenID") %>% 
               mutate(probenID = as.numeric(probenID)))


## T1 protein plot -------------------
selectedPro <- "OID20722"
abTiter <- "H1N1_titerFC"

plot_dat <- inputDat %>% filter(time == "T1") %>% 
  select(patientID, condition, disease, category, {{abTiter}},
         OID20481, OID20624, OID20765) # IL17D, TNFSF12, CCL22
cowplot::plot_grid(
  ggplot(plot_dat, aes(x = H1N1_titerFC, y = OID20481)) +
    geom_point(aes(color = category, shape = condition), 
               size = 3, alpha = 0.8,
               position = position_jitter(w = 0.08, h = 0.1)) +
    geom_smooth(method = "lm", se = FALSE) + theme_bw(),
  ggplot(plot_dat, aes(x = H1N1_titerFC, y = OID20624)) +
    geom_point(aes(color = category, shape = condition), 
               size = 3, alpha = 0.8,
               position = position_jitter(w = 0.08, h = 0.1)) +
    geom_smooth(method = "lm", se = FALSE) + theme_bw(),
  ggplot(plot_dat, aes(x = H1N1_titerFC, y = OID20765)) +
    geom_point(aes(color = category, shape = condition), 
               size = 3, alpha = 0.8,
               position = position_jitter(w = 0.08, h = 0.1)) +
    geom_smooth(method = "lm", se = FALSE) + theme_bw(),
  nrow = 3
)

plot_dat$category <- factor(plot_dat$category, levels = c("NR", "Other", "TR"))
plot_dat$disease <- factor(plot_dat$disease, 
                           levels = c("healthy", "cirrhosis"))
plot_dat$condition <- factor(plot_dat$condition, 
                             levels = c("healthy", "compensated cirrhosis", "decompensated cirrhosis"))
cowplot::plot_grid(
  ggboxplot(plot_dat, x = "category", y = "OID20481", palette = "jco", add = "jitter"),
  ggboxplot(plot_dat, x = "disease", y = "OID20481", palette = "jco", add = "jitter"),
  ggboxplot(plot_dat, x = "condition", y = "OID20481", palette = "jco", add = "jitter") + 
    rotate_x_text(20),
  ggboxplot(plot_dat, x = "category", y = "OID20624", palette = "jco", add = "jitter"),
  ggboxplot(plot_dat, x = "disease", y = "OID20624", palette = "jco", add = "jitter"),
  ggboxplot(plot_dat, x = "condition", y = "OID20624", palette = "jco", add = "jitter") + 
    rotate_x_text(20),
  ggboxplot(plot_dat, x = "category", y = "OID20765", palette = "jco", add = "jitter"),
  ggboxplot(plot_dat, x = "disease", y = "OID20765", palette = "jco", add = "jitter"),
  ggboxplot(plot_dat, x = "condition", y = "OID20765", palette = "jco", add = "jitter") + 
    rotate_x_text(20),
  nrow = 3, byrow = FALSE
)






