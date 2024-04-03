library(readxl)
library(tidyverse)

# cirrhosis donor information ------------------------------
cirrhosis_rawDat <- read_xlsx("../metadata/Cohort information Zirrflu Nhan.xlsx",
               sheet = "cirrhotic patients", col_names = FALSE)

cirrhosis_temp <- cirrhosis_rawDat %>% 
  mutate(season = c("NA", "season", rep("2019", 32), "NA", "NA", rep("2020", 9))) %>%
  select(-c(1:5, 17:18)) %>% relocate(season)

colnames(cirrhosis_temp) <- cirrhosis_temp[2, ]

cirrhosis_Info <- cirrhosis_temp %>% 
  slice(-c(1, 2, 35, 36))

# healthy control information ---------------------------------
control_rawDat <- read_xlsx("../metadata/Cohort information Zirrflu Nhan.xlsx",
               sheet = "healthy controls", col_names = FALSE)

control_temp <- control_rawDat %>% slice(-2) %>%
  fill(1, 4) %>% select(-c(2, 3))

colnames(control_temp) <- c("season", "center", control_temp[1,-c(1:2)])

control_Info <- control_temp %>% slice(-1)
control_Info$`previous flu vaccination`[59:69] <- rep("please ask Janyn for the selected patients")

# merge data together ------------------------------------------
all_info <- cirrhosis_Info %>% 
  full_join(control_Info, 
            by = c("season", "Zirrflu-number" = "Zirrflu number" ,
                   "age (years)", "height (cm)", "height (m)", 
                   "weight (kg)", "BMI", "previous flu vaccination"))


# check the correlation of cirrhosis state with vaccine response --------------
a <- ZirFlu$HAItiter %>% filter(disease == "cirrhosis") %>%
  left_join(cirrhosis_Info, by = c("patientID" = "Zirrflu-number", "season")) %>%
  mutate(BMI = as.numeric(BMI),
         `CHILD-Pugh ponts` = as.numeric(`CHILD-Pugh ponts`),
         `CHILD-Pugh state of cirrhosis (A=1, B=2, C=3)` = as.numeric(`CHILD-Pugh state of cirrhosis (A=1, B=2, C=3)`),
         `MELD points` = as.numeric(`MELD points`))

cor(a$H1N1_abFC, as.numeric(a$BMI), use = "complete.obs")
cor(a$H1N1_abFC, as.numeric(a$`CHILD-Pugh ponts`), use = "complete.obs")
cor(a$H1N1_abFC, as.numeric(a$`CHILD-Pugh state of cirrhosis (A=1, B=2, C=3)`), use = "complete.obs")
cor(a$H1N1_abFC, as.numeric(a$`MELD points`), use = "complete.obs")

cor(a$H3N2_abFC, as.numeric(a$`CHILD-Pugh ponts`), use = "complete.obs")
cor(a$H3N2_abFC, as.numeric(a$`CHILD-Pugh state of cirrhosis (A=1, B=2, C=3)`), use = "complete.obs")
cor(a$H3N2_abFC, as.numeric(a$`MELD points`), use = "complete.obs")

ggplot(a, aes(x = H1N1_abFC, y = `CHILD-Pugh ponts`)) +
  geom_point(position = position_jitter(w = 0.08, h = 0.1)) +
  geom_smooth(method = "lm", se = FALSE) + theme_bw()

library(ggpubr)
a %>% 
  ggscatter(x = "H1N1_abFC", y = "CHILD-Pugh ponts", 
            color = "condition", add = "reg.line") +
  stat_cor(aes(color = condition), method = "pearson")

a %>% 
  ggscatter(x = "H3N2_abFC", y = "CHILD-Pugh ponts", 
            color = "condition", add = "reg.line") +
  stat_cor(aes(color = condition), method = "pearson")

a %>% 
  ggscatter(x = "Byamagata_abFC", y = "CHILD-Pugh ponts", 
            color = "condition", add = "reg.line") +
  stat_cor(aes(color = condition), method = "pearson")

a %>% 
  ggscatter(x = "Bvictoria_abFC", y = "CHILD-Pugh ponts", 
            color = "condition", add = "reg.line") +
  stat_cor(aes(color = condition), method = "pearson")

a$esophageal_varices <- c(0, NA, 1, 1, NA, 1, 1, 0, 1, 1, 1, 0, 
                          1, 1, 0, NA, 1, 1, 1, 1, 1, 1, 1, 1, 
                          0, 1, 1, 1, NA, 1, 0, 1, 1, 1, 1, 1)

variable <- "CHILD-Pugh ponts"
variable <- "CHILD-Pugh state of cirrhosis (A=1, B=2, C=3)"
variable <- "MELD points"
variable <- "esophageal_varices"
cowplot::plot_grid(
  a %>% 
  ggscatter(x = "H1N1_abFC", y = variable, 
            color = "condition", add = "reg.line") +
  stat_cor(aes(color = condition), method = "pearson"),
  a %>% 
  ggscatter(x = "H3N2_abFC", y = variable, 
            color = "condition", add = "reg.line") +
  stat_cor(aes(color = condition), method = "pearson"),
  a %>% 
  ggscatter(x = "Byamagata_abFC", y = variable, 
            color = "condition", add = "reg.line") +
  stat_cor(aes(color = condition), method = "pearson"),
  a %>% 
  ggscatter(x = "Bvictoria_abFC", y = variable, 
            color = "condition", add = "reg.line") +
  stat_cor(aes(color = condition), method = "pearson"))



