all_res <- readRDS("../../output/gridtests.rds")

all_res <- all_res %>%
  mutate(llik_diff = pmax(diff, 0),
         n = factor(n),
         family = factor(family, levels = c("normalmix", "unimix", "npmle"))) %>%
  mutate(family = recode_factor(
    family, 
    normalmix = "Scale mixtures of normals",
    unimix = "Symmetric unimodal priors",
    npmle = "NPMLE"
  ))

plt <- ggplot(all_res, aes(x = n, y = llik_diff)) +
  geom_boxplot() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  stat_summary(fun = mean, geom = "point", size = 3, shape = 8) +
  theme_minimal() +
  labs(x = "\nNumber of observations",
       y = "Difference in log likelihood from optimal\n") +
  facet_wrap(~family) +
  theme(strip.text = element_text(size = 12))
  
ggsave(paste0("../../figs/gridtest_all.png"), height = 5, width = 7.5)


# for (ebnm_family in unique(all_res$family)) {
#   plt <- ggplot(all_res %>% filter(family == ebnm_family), aes(x = n, y = llik_diff)) +
#     geom_boxplot() +
#     geom_hline(yintercept = 1, linetype = "dashed") +
#     theme_minimal() +
#     labs(x = "\nNumber of observations",
#          y = "Difference in log likelihood from optimal\n")
#   
#   ggsave(paste0("../../figs/gridtest_", ebnm_family, ".png"), height = 3.5, width = 5)
# }