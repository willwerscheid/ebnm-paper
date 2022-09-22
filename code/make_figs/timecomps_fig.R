## Timing different prior families.

res <- readRDS("../../output/timecomps.rds")

ggplot(res, aes(x = n, y = t, col = fn)) +
  geom_point() +
  geom_line(size = 0.2) +
  scale_x_log10() +
  scale_y_log10() +
  labs(x = "\nNumber of observations", y = "Time elapsed (s)\n", col = "ebnm function") +
  theme_minimal()

ggsave("../../figs/timecomps.png", width = 6, height = 3.375)
