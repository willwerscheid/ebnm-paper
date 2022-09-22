## Performance comparisons: ebnm_npmle vs. REBayes

rebayes <- readRDS("../../output/rebayes.rds")
rebayes <- rebayes %>%
  mutate(homosked = ifelse(homosked, "homoskedastic", "heteroskedastic"),
         homosked = fct_relevel(homosked, "homoskedastic", "heteroskedastic")) %>%
  group_by(n, n_gridpts, homosked) %>%
  mutate(t_penalty = log2(mean_t / min(mean_t))) %>%
  ungroup()

ggplot(rebayes, aes(x = n, y = package, fill = t_penalty)) +
  scale_x_log10() +
  scale_y_discrete(limits = c("REBayes", "ebnm")) +
  geom_raster() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 1, 
                       limits = c(0, 4)) +
  labs(x = "\nn", y = "package\n", fill = "log2(time/\nfastest time)\n") +
  facet_grid(cols = vars(n_gridpts), rows = vars(homosked))
  
ggsave("../../figs/rebayes.png", width = 7, height = 5)
