library(tidyverse)
library(rstatix)
library(ggpubr)


# overview ------------------------
dat <- ZirFlu$HAItiter_2019_2020 %>% mutate(season = "2019") %>%
  full_join(ZirFlu$HAItiter_2020_2021 %>% mutate(season = "2020"))

dat %>% count(season)
dat %>% filter(season == "2019") %>% count(condition)
dat %>% filter(season == "2020") %>% count(condition)

# NOTE:
# For the ZirFlu cohort 2019-2020, we excluded the participant Z-62, 
# because unfortunately no serum samples were available in the biobank.

# The participant with the ID in 2020/2021: Z-01-99-069 has azathioprine in his medication 
# and should be excluded from the analysis as well. 


# box plot -------------

dat1 <- ZirFlu$HAItiter_2019_2020 %>% 
  select(patientID, matches("ab")) %>%
  pivot_longer(!patientID, names_to = "strains", values_to = "ab_value") %>%
  mutate(season = "2019-2020")

dat2 <- ZirFlu$HAItiter_2020_2021 %>% 
  select(patientID, matches("ab")) %>%
  pivot_longer(!patientID, names_to = "strains", values_to = "ab_value") %>%
  mutate(season = "2020-2021")

dat3 <- full_join(dat1, dat2) %>% as.data.frame() %>%
  mutate(ab_value = as.numeric(ab_value)) %>%
  mutate(strains = as.factor(strains), season = as.factor(season))

dat4 <- dat3 %>% full_join(ZirFlu$donorInfo2)
dat4$condition <- factor(dat4$condition, 
                        levels = c("healthy", "compensated cirrhosis", "decompensated cirrhosis" ))

dat4 %>%
  ggboxplot(x = "condition", y = "ab_value",
            color = "black", fill = "season", palette = "ucscgb", 
            bxp.errorbar = TRUE, 
            outlier.shape = 1, outlier.size = 4,outlier.color = "grey")


bxplot <- dat4 %>%
  ggboxplot(x = "condition", y = "ab_value",
            color = "black", fill = "season", palette = "npg", 
            bxp.errorbar = TRUE, 
            outlier.shape = 1, outlier.size = 4, outlier.color = "grey")

stat.test <- dat4 %>%
  group_by(condition) %>% t_test(ab_value ~ season) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj") %>% 
  add_xy_position(x = "condition", dodge = 0.08)

bxplot + 
  stat_pvalue_manual(stat.test, label = "p", tip.length = 0)

stat.test2 <- dat4 %>% 
  t_test(ab_value ~ condition, p.adjust.method = "bonferroni") %>%
  add_xy_position(x = "condition")
stat.test2

bxplot + 
  #stat_pvalue_manual(stat.test, label = "p", tip.length = 0.01) +
  stat_pvalue_manual(stat.test2, label = "p", tip.length = 0.02, step.increase = 0.05)

# option 2
bxp <- dat4 %>%
  ggboxplot(x = "season", y = "ab_value", fill = "condition",
            bxp.errorbar = TRUE, palette = "npg",
            outlier.shape = 1, outlier.size = 4, outlier.color = "grey")


stat.test <- dat4 %>%
  group_by(season) %>%
  t_test(ab_value ~ condition)%>%
  add_xy_position(x = "season", dodge = 0.8)
bxp + 
  stat_pvalue_manual(
    stat.test, label = "p.adj", tip.length = 0.01,
    bracket.nudge.y = 1
  ) # + scale_y_continuous(expand = expansion(mult = c(0, 0.15)))

# group line plots -------------
rm(dat1, dat2, dat3, dat4)

dat1 <- ZirFlu$HAItiter_2019_2020 %>% 
  select(-matches("ab"), -c(vaccine_response, category, condition)) %>%
  pivot_longer(!patientID, names_to = "strains_day", values_to = "HAI_value") %>%
  mutate(season = "2019-2020") %>%
  separate(strains_day, c("strains", "day")) # need to by fix the date split

dat2 <- ZirFlu$HAItiter_2020_2021 %>% 
  select(-matches("ab"), -c(vaccine_response, category, condition)) %>%
  pivot_longer(!patientID, names_to = "strains_day", values_to = "HAI_value") %>%
  mutate(season = "2020-2021") %>%
  separate(strains_day, c("strains", "day")) # need to by fix the date split

dat3 <- full_join(dat1, dat2) %>% as.data.frame() %>%
  mutate(day = ifelse(day == "d0", "T1", ifelse(day == "d21", "T2", "T3")))

dat4 <- dat3 %>% full_join(ZirFlu$donorInfo2)
dat4$condition <- factor(dat4$condition, 
                         levels = c("healthy", "compensated cirrhosis", "decompensated cirrhosis" ))
datplot <- dat4 %>% filter(season == "2019-2020")
datplot <- dat4 %>% filter(season == "2019-2020") %>%
  mutate(condition = ifelse(condition == "healthy", "healthy", "cirrhosis"))

lplot <- datplot %>% ggline(x = "day", y = "HAI_value", add = "mean_sd", 
                color = "condition", palatte = "npg")

stat.test <- datplot %>%
  group_by(condition) %>% t_test(HAI_value ~ day) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj") %>% 
  add_xy_position(x = "condition", fun = "mean_sd")
lplot + stat_pvalue_manual(stat.test, label = "p.adj.signif",
                           linetype = "blank")

pwc <- datplot %>%
  group_by(condition) %>%
  t_test(HAI_value ~ day, p.adjust.method = "bonferroni") %>%
  add_xy_position(x = "day", fun = "mean_sd", dodge = 0.8)
lplot + stat_pvalue_manual(pwc, color = "condition", step.group.by = "condition",
                          tip.length = 0, step.increase = 0.1,
                          bracket.nudge.y = 2
                          )  +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)))
# HAi repsonse to age -----------
rm(dat1, dat2, dat3, dat4)

dat1 <- ZirFlu$HAItiter_2019_2020 %>% 
  select(patientID, matches("ab")) %>%
  pivot_longer(!patientID, names_to = "strains", values_to = "ab_value") %>%
  mutate(season = "2019-2020")

dat2 <- ZirFlu$HAItiter_2019_2020 %>% 
  full_join(ZirFlu$donorInfo, by = c("patientID", "condition"))

dat2 %>% 
  ggscatter(x = "age", y = "H1N1_abFC", 
            add = "reg.line") +
  stat_cor(method = "pearson", label.x = 60, label.y = 7)

dat2 %>% 
  ggscatter(x = "age", y = "H1N1_abFC", color = "condition",
            add = "reg.line") +
  stat_cor(aes(color = condition), method = "pearson", label.x = 60)
# calculate the correlation between A/H1N1, A/H3N2, B/viectoria, b/Yamagata
cor(ZirFlu$HAItiter_2019_2020$H1N1_abFC, ZirFlu$HAItiter_2019_2020$H3N2_abFC)
cowplot::plot_grid(
  ZirFlu$HAItiter_2019_2020 %>%
  ggscatter(x = "H1N1_abFC", y = "H3N2_abFC",
            add = "reg.line") +
  stat_cor(method = "pearson", label.x = 6, label.y = 7.5),
  ZirFlu$HAItiter_2019_2020 %>%
  ggscatter(x = "Bvictoria_abFC", y = "Byamagata_abFC",
            add = "reg.line") +
  stat_cor(method = "pearson", label.x = 6, label.y = 7.5),
)
ZirFlu$HAItiter_2019_2020 %>%
  ggscatter(x = "H1N1_abFC", y = "Bvictoria_abFC",
            add = "reg.line") +
  stat_cor(method = "pearson", label.x = 6, label.y = 7)
ZirFlu$HAItiter_2019_2020 %>%
  ggscatter(x = "H3N2_abFC", y = "Bvictoria_abFC",
            add = "reg.line") +
  stat_cor(method = "pearson", label.x = 6, label.y = 7)
