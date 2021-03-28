## Performance comparisons: ebnm_normal_scale_mixture vs. ashr.

ashr <- readRDS("../../output/ebnmvashr.rds")
ashr <- ashr %>%
  mutate(prior = case_when(
    sim_fn == "null_sim" ~ "prior_null",
    sim_fn == "pt_sim" ~ "prior_true",
    TRUE ~ "prior_misspec"
  )) %>%
  mutate(homosked = ifelse(homosked, "homoskedastic", "heteroskedastic"),
         homosked = fct_relevel(homosked, "homoskedastic", "heteroskedastic"),
         prior = fct_relevel(prior, "prior_true", "prior_null", "prior_misspec")) %>%
  group_by(n, homosked, prior) %>%
  mutate(t_penalty = log2(mean_t / min(mean_t))) %>%
  ungroup()

ggplot(ashr, aes(x = n, y = package, fill = t_penalty)) +
  scale_x_log10() +
  scale_y_discrete(limits = c("ashr", "ebnm")) +
  geom_raster() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 1, 
                       limits = c(0, 4)) +
  labs(x = "\nn", y = "package\n", fill = "t vs. fastest\n(log2 scale)\n") +
  facet_grid(cols = vars(prior), rows = vars(homosked))
  
ggsave(paste0("../../figs/ebnmvashr.png"), width = 6, height = 3.375)
