all_res <- readRDS("../../output/gridtests.rds")

all_res <- all_res %>%
  mutate(llik_diff = pmax(diff, 0),
         n = factor(n))

for (ebnm_family in unique(all_res$family)) {
  plt <- ggplot(all_res %>% filter(family == ebnm_family), aes(x = n, y = llik_diff)) +
    geom_boxplot() +
    geom_hline(yintercept = 1, linetype = "dashed") +
    theme_minimal() +
    labs(x = "Number of observations",
         y = "Difference in log likelihood from optimal")
  
  ggsave(paste0("../../figs/gridtest_", ebnm_family, ".png"), height = 3.5, width = 5)
}
