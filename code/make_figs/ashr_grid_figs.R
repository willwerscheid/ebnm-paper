# Quality of ashr grid approximations.

tib <- readRDS("../../output/smngrid.rds")
ggplot(tib, aes(x = m, y = ub)) +
  geom_line() +
  scale_y_log10() +
  labs(x = "\ngrid multiplier", y = "KL divergence\n") +
  theme_minimal()
ggsave("../../figs/smnKLdiv.png", width = 6, height = 3.375)

tib <- readRDS("../../output/symmgrid.rds")
ggplot(tib, aes(x = min_space, y = KL)) + 
  geom_line() +
  scale_y_log10() +
  scale_x_continuous(labels = function(x) paste0(x, "s")) + 
  labs(x = "\ngrid spacing", y = "KL-divergence\n") +
  theme_minimal()
ggsave("../../figs/symmKLdiv.png", width = 6, height = 3.375)
