set.seed(666)

n <- 8
x <- c(28, 8, -3, 7, -1, 1, 18, 12)
s <- c(15, 10, 16, 11, 9, 11, 10, 18)

ebnm_res_mode0 <- ebnm::ebnm_point_normal(
  x = x, s = s, mode = 0, output = output_all()
)
samp <- ebnm_res_mode0$posterior_sampler(10000)
CI <- apply(samp, 2, quantile, c(.025, .975))

tib <- tibble(
  idx = 1:n,
  x = x, 
  s = s, 
  theta_hat = ebnm_res_mode0$posterior$mean,
  CI_lower = CI[1, ],
  CI_upper = CI[2, ]
)

plt <- ggplot(tib, aes(x = idx, y = x)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = x - 2 * s, ymax = x + 2 * s)) +
  geom_point(aes(x = idx, y = theta_hat, col = "red"), size = 3) +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper, col = "red")) +
  theme_minimal() +
  labs(x = "", y = "Treatment effect") +
  theme(legend.position = "none",
        axis.text.x = element_blank())

ggsave("../../figs/eight_schools_mode0.png", height = 3.5, width = 5)

ebnm_res_estmode <- ebnm::ebnm_point_normal(
  x = x, s = s, mode = "estimate", output = output_all()
)
samp <- ebnm_res_estmode$posterior_sampler(10000)
CI <- apply(samp, 2, quantile, c(.025, .975))

tib <- tibble(
  idx = 1:n,
  x = x, 
  s = s, 
  theta_hat = ebnm_res_estmode$posterior$mean,
  CI_lower = CI[1, ],
  CI_upper = CI[2, ]
)

plt <- ggplot(tib, aes(x = idx, y = x)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = x - 2 * s, ymax = x + 2 * s)) +
  geom_point(aes(x = idx, y = theta_hat, col = "red"), size = 3) +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper, col = "red")) +
  theme_minimal() +
  labs(x = "", y = "Treatment effect") +
  theme(legend.position = "none",
        axis.text.x = element_blank())

ggsave("../../figs/eight_schools_modeest.png", height = 3.5, width = 5)
