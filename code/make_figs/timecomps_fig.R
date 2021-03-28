## Timing different prior families.

res <- readRDS("../../output/timecomps.rds")

ggplot(res, aes(x = n, y = t, col = fn)) +
  geom_point() +
  scale_x_log10() +
  scale_y_log10() +
  labs(x = "number of observations", y = "time elapsed (s)", col = "ebnm function")

ggsave("../../figs/timecomps.png", width = 6, height = 3.375)
