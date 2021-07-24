# Quality of ashr grid approximations.

tib <- ebnm:::smngrid
ggplot(tib, aes(x = m, y = ub)) +
  geom_line() +
  scale_y_log10() +
  labs(x = "\ngrid multiplier", y = "KL divergence\n") +
  theme_minimal()
ggsave("../../figs/smnKLdiv.png", width = 6, height = 3.375)
