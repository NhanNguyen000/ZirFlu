# ggplot: https://www.datanovia.com/en/lessons/how-to-do-a-t-test-in-r-calculation-and-reporting/how-to-do-paired-t-test-in-r/
ZirFlu$protein_annot[which(ZirFlu$protein_annot$OlinkID == "OID20507"),]
library(ggpubr)
ggboxplot(abTiter_ttest, x = "group", y = "H1N1_d21.35") + 
  stat_pvalue_manual(meta_ttest_outcome$abTiter[2,] %>% add_xy_position(x = "group"), 
                     tip.length = 0) +
  labs(subtitle = get_test_label(meta_ttest_outcome$abTiter[2,], detailed = TRUE))

ggpaired(cytokine_ttest, x = "Measure", y = "GRO....CXCL1.", # have issue with its
         order = c("d0", "d21-35")) + 
  stat_pvalue_manual(meta_ttest_outcome$cytokine_pairedTtest[2,] %>% add_xy_position(x = "group"), 
                     tip.length = 0) +
  labs(subtitle = get_test_label(meta_ttest_outcome$cytokine_pairedTtest[2,], detailed = TRUE))
