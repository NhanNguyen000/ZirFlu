library(tidyverse)
library(rstatix)
library(ggpubr)


# overview ---------------------------------------------------------------------
ZirFlu$HAItiter %>% count(season)
ZirFlu$HAItiter %>% filter(season == "2019") %>% count(condition)
ZirFlu$HAItiter %>% filter(season == "2020") %>% count(condition)

# NOTE:
# The participant Z-62 season 2019 were excluded because no serum samples were available.
# The participant Z-01-99-069 season 2020 has azathioprine in his medication 
# and were excluded from the analysis.

# box plot (abFC)----------------------------------------------------------------------
HAI_ab <- ZirFlu$HAItiter %>% select(patientID, season, condition,  matches("ab")) %>%
  pivot_longer(!c("patientID", "season", "condition"), 
               names_to = "strains", values_to = "ab_value") %>%
  mutate(condition = factor(condition, 
                            levels = c("healthy", "compensated cirrhosis", "decompensated cirrhosis"))) %>%
  mutate(disease = ifelse(condition == "healthy", "healthy", "cirrhosis")) %>%
  mutate(disease = factor(disease,
                          levels = c("healthy", "cirrhosis")))

# option 1.1: group per season, healthy vs. 2 cirrhosis
bxp <- HAI_ab %>%
  ggboxplot(x = "season", y = "ab_value", fill = "condition",
            bxp.errorbar = TRUE, palette = "npg",
            outlier.shape = 1, outlier.size = 4, outlier.color = "grey")

stat.test <- HAI_ab %>% group_by(season) %>%
  t_test(ab_value ~ condition)%>% add_xy_position(x = "season", dodge = 0.8)

bxp + 
  stat_pvalue_manual(
    stat.test, label = "p.adj", tip.length = 0.01, bracket.nudge.y = 1)  + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)))

# option 1.2: group per season, healthy vs. 1 cirrhosis
bxp <- HAI_ab %>%
  ggboxplot(x = "season", y = "ab_value", fill = "disease",
            bxp.errorbar = TRUE, palette = "npg",
            outlier.shape = 1, outlier.size = 4, outlier.color = "grey")

stat.test <- HAI_ab %>% group_by(season) %>%
  t_test(ab_value ~ disease)%>% add_xy_position(x = "season", dodge = 0.8)

bxp + 
  stat_pvalue_manual(
    stat.test, label = "p", tip.length = 0.01, bracket.nudge.y = 1)  + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)))

# option 2:  group per healthy condition, season 2019 vs. 2020
HAI_ab %>%
  ggboxplot(x = "condition", y = "ab_value",
            color = "black", fill = "season", palette = "ucscgb", 
            bxp.errorbar = TRUE, 
            outlier.shape = 1, outlier.size = 4,outlier.color = "grey")

bxplot <- HAI_ab %>%
  ggboxplot(x = "condition", y = "ab_value",
            color = "black", fill = "season", palette = "npg", 
            bxp.errorbar = TRUE, 
            outlier.shape = 1, outlier.size = 4, outlier.color = "grey")

stat.test <- HAI_ab %>%
  group_by(condition) %>% t_test(ab_value ~ season) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj") %>% 
  add_xy_position(x = "condition", dodge = 0.08)

bxplot + 
  stat_pvalue_manual(stat.test, label = "p", tip.length = 0)

stat.test2 <- HAI_ab %>% 
  t_test(ab_value ~ condition, p.adjust.method = "bonferroni") %>%
  add_xy_position(x = "condition")

bxplot + 
  #stat_pvalue_manual(stat.test, label = "p", tip.length = 0.01) +
  stat_pvalue_manual(stat.test2, label = "p", tip.length = 0.02, step.increase = 0.05)


# group line plots -----------------------------------------------------------------
HAI_level <- ZirFlu$HAItiter %>% select(-matches("ab"), -c(vaccine_response, category)) %>%
  pivot_longer(!c("patientID", "condition", "season"), 
               names_to = "strains_day", values_to = "HAI_value") %>%
  separate(strains_day, c("srtains", "day")) %>%
  mutate(disease = ifelse(condition == "healthy", "healthy", "cirrhosis"))

datplot <- HAI_level %>% filter(season == "2019")
datplot <- HAI_level %>% filter(season == "2020")
datplot <- HAI_level
# option 1.1: group per 3 healthy conditions, compare per time point
lineplot <-  datplot %>% 
  ggline(x = "day", y = "HAI_value", add = "mean_sd", 
         color = "condition", palatte = "npg")

stat.test <- datplot %>%
  group_by(condition) %>% t_test(HAI_value ~ day) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj") %>% 
  add_xy_position(x = "condition", fun = "mean_sd")

lineplot + stat_pvalue_manual(stat.test, label = "p.adj.signif",
                           linetype = "blank")

# option 1.2: group per 2 healthy conditions, compare between time point
lineplot <-  datplot %>% 
  ggline(x = "day", y = "HAI_value", add = "mean_sd", 
         color = "disease", palatte = "npg")

stat.test <- datplot %>%
  group_by(disease) %>% t_test(HAI_value ~ day) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj") %>% 
  add_xy_position(x = "disease", fun = "mean_sd")

pwc <- datplot %>%
  group_by(disease) %>%
  t_test(HAI_value ~ day, p.adjust.method = "bonferroni") %>%
  add_xy_position(x = "day", fun = "mean_sd", dodge = 0.8)

lineplot + stat_pvalue_manual(pwc, color = "disease", step.group.by = "disease",
                           tip.length = 0, step.increase = 0.1,
                           bracket.nudge.y = 2)  +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)))

# calculate the correlation between A/H1N1, A/H3N2, B/Victoria, b/Yamagata ------------
library(corrplot)

# 16 identical patients in both 2 season:
patient_in2seasons <- ZirFlu$HAItiter_2019 %>% 
  select(patientID, condition,  matches("ab")) %>% rename_with(~sub("abFC", "abFC_2019", .x)) %>%
  inner_join(ZirFlu$HAItiter_2020 %>% 
               select(patientID, condition,  matches("ab"))  %>% rename_with(~sub("abFC", "abFC_2020", .x)))

corDat_2seasons <- patient_in2seasons %>% select(-condition) %>% 
  column_to_rownames("patientID") %>% as.matrix() %>% cor()
testRes = cor.mtest(corDat_2seasons, conf.level = 0.95)

# corrplot 
corDat_2seasons %>% corrplot()
corDat_2seasons %>% corrplot(type = "upper", tl.col = "black", 
                             tl.srt = 45, addCoef.col = "black")

p2 <- corDat_2seasons %>% corrplot(type = "upper", 
                             tl.col = "black", tl.srt = 45,
                             addCoef.col = "black", insig = "blank",
                             p.mat = testRes$p, sig.level = 0.05)
text(p2$corrPos$x, p2$corrPos$y, round(p2$corrPos$corr, 2))

# per season
dat <- ZirFlu$HAItiter_2019
dat <- ZirFlu$HAItiter_2020

corDat <- dat %>% select(matches("ab")) %>%  as.matrix() %>% cor()
testRes = cor.mtest(corDat, conf.level = 0.95)

corDat %>% corrplot()
corDat %>% corrplot(type = "upper", tl.col = "black", 
                   tl.srt = 45, addCoef.col = "black")

p2 <- corDat %>% corrplot(type = "upper", tl.col = "black", tl.srt = 60,
                          addCoef.col = "black", insig = "blank",
                          p.mat = testRes$p, sig.level = 0.05)
text(p2$corrPos$x, p2$corrPos$y, round(p2$corrPos$corr, 2))

# other

cor(ZirFlu$HAItiter_2019$H1N1_abFC, ZirFlu$HAItiter_2019$H3N2_abFC)
cowplot::plot_grid(
  ZirFlu$HAItiter_2019 %>%
    ggscatter(x = "H1N1_abFC", y = "H3N2_abFC",
              add = "reg.line") +
    stat_cor(method = "pearson", label.x = 6, label.y = 7.5),
  ZirFlu$HAItiter_2019 %>%
    ggscatter(x = "Bvictoria_abFC", y = "Byamagata_abFC",
              add = "reg.line") +
    stat_cor(method = "pearson", label.x = 6, label.y = 7.5),
)
ZirFlu$HAItiter_2019 %>%
  ggscatter(x = "H1N1_abFC", y = "Bvictoria_abFC",
            add = "reg.line") +
  stat_cor(method = "pearson", label.x = 6, label.y = 7)
ZirFlu$HAItiter_2019 %>%
  ggscatter(x = "H3N2_abFC", y = "Bvictoria_abFC",
            add = "reg.line") +
  stat_cor(method = "pearson", label.x = 6, label.y = 7)

# HAI heterogenitety - 4 strains (need to check more)-----------------------------------
dat <- HAI_ab
ggplot(data=HAI_ab, aes(y=patientID, x = ab_value, group = strains)) +
  geom_line()+
  geom_point()

# HAi repsonse to age -------------------------------------------------------
rm(dat1, dat2, dat3, dat4)

dat1 <- ZirFlu$HAItiter_2019 %>% 
  select(patientID, matches("ab")) %>%
  pivot_longer(!patientID, names_to = "strains", values_to = "ab_value") %>%
  mutate(season = "2019-2020")

dat2 <- ZirFlu$HAItiter_2019 %>% 
  full_join(ZirFlu$donorInfo, by = c("patientID", "condition"))

# H1N1 
dat2 %>% 
  ggscatter(x = "age", y = "H1N1_abFC", add = "reg.line") +
  stat_cor(method = "pearson", label.x = 60, label.y = 7)

dat2 %>% 
  ggscatter(x = "age", y = "H1N1_abFC", color = "condition", add = "reg.line") +
  stat_cor(aes(color = condition), method = "pearson", label.x = 60)

dat2 %>% 
  ggscatter(x = "age", y = "H1N1_abFC", color = "disease", add = "reg.line") +
  stat_cor(aes(color = disease), method = "pearson", label.x = 60)

# H3N2
strain_abFC <- "H3N2_abFC"
strain_abFC <- "Byamagata_abFC"
strain_abFC <- "Bvictoria_abFC"
cowplot::plot_grid(
  dat2 %>% 
    ggscatter(x = "age", y = strain_abFC, add = "reg.line") +
    stat_cor(method = "pearson", label.x = 60, label.y = 7.5),
  dat2 %>% 
    ggscatter(x = "age", y = strain_abFC, color = "disease", add = "reg.line") +
    stat_cor(aes(color = disease), method = "pearson", label.x = 60)
)

