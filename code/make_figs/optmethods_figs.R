## Performance comparisons: ebnm optmethods.

optmethods <- readRDS("../../output/optmethods.rds")
optmethods <- optmethods %>%
  mutate(prior = case_when(
    sim_fn == "null_sim" ~ "prior_null",
    sim_fn %in% c("pn_sim", "pl_sim") ~ "prior_true",
    ebnm_fn %in% c("ebnm_pn_estmode", "ebnm_pl_estmode") &
      sim_fn %in% c("pn_nzmode_sim", "pl_nzmode_sim") ~ "prior_true",
    TRUE ~ "prior_misspec"
  )) %>%
  mutate(homosked = ifelse(homosked, "homoskedastic", "heteroskedastic"),
         homosked = fct_relevel(homosked, "homoskedastic", "heteroskedastic"),
         prior = fct_relevel(prior, "prior_true", "prior_null", "prior_misspec")) %>%
  group_by(ebnm_fn, n, homosked, prior) %>%
  mutate(t_penalty = log2(mean_t / min(mean_t))) %>%
  ungroup()

for (next_ebnm_fn in unique(optmethods$ebnm_fn)) {
  plt_title <- ifelse(next_ebnm_fn %in% c("ebnm_point_normal", "ebnm_point_laplace"),
                      "Mode: fixed at zero\n",
                      "Mode: estimated\n")
  
  plt <- ggplot(optmethods %>% filter(ebnm_fn == next_ebnm_fn), 
                aes(x = n, y = optmethod, fill = t_penalty)) +
    scale_x_log10() +
    scale_y_discrete(limits = c("nograd_lbfgsb", "nograd_nlm",
                                "lbfgsb", "nohess_nlm",
                                "trust", "nlm")) +
    geom_raster() +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 1, 
                         limits = c(0, 4)) +
    labs(x = "\nn", y = "optmethod\n", fill = "log2(time/\nfastest time)\n",
         title = plt_title) +
    facet_grid(cols = vars(prior), rows = vars(homosked))

  ggsave(paste0("../../figs/optmethods_", next_ebnm_fn, ".png"), width = 7, height = 5,
         plot = plt)
}
