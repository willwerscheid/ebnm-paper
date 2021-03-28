tib <- readRDS("../../output/optgrids.rds")
ggplot(tib, aes(x = grid, y = KL)) + 
  geom_point() +
  geom_text(aes(label = ifelse(idx %% 5 == 0, paste0(idx, "th"), '')),
            hjust = 0,
            vjust = 2,
            size = 2) +
  scale_y_log10() +
  scale_x_continuous(labels = function(x) paste0(x, "s")) + 
  labs(x = "\ngrid points", y = "KL-divergence\n") +
  theme_minimal()
ggsave("../../figs/optgrids.png", width = 6, height = 3.375)
