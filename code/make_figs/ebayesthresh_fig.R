## Performance comparisons: ebnm_point_laplace vs. EBayesThresh.

ebayesthresh <- readRDS("../../output/ebayesthresh.rds")
ebayesthresh <- ebayesthresh %>%
  mutate(prior = case_when(
    sim_fn == "null_sim" ~ "prior_null",
    sim_fn == "pl_sim" ~ "prior_true",
    TRUE ~ "prior_misspec"
  )) %>%
  mutate(homosked = ifelse(homosked, "homoskedastic", "heteroskedastic"),
         homosked = fct_relevel(homosked, "homoskedastic", "heteroskedastic"),
         prior = fct_relevel(prior, "prior_true", "prior_null", "prior_misspec")) %>%
  group_by(n, homosked, prior) %>%
  mutate(t_penalty = log2(mean_t / min(mean_t))) %>%
  ungroup()

ggplot(ebayesthresh, aes(x = n, y = package, fill = t_penalty)) +
  scale_x_log10() +
  scale_y_discrete(limits = c("EBayesThresh", "ebnm")) +
  geom_raster() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 1, 
                       limits = c(0, 4)) +
  labs(x = "\nn", y = "package\n", fill = "log2(time/\nfastest time)\n") +
  facet_grid(cols = vars(prior), rows = vars(homosked))
  
ggsave("../../figs/ebayesthresh.png", width = 7, height = 5)
