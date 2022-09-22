library(tidyverse)
library(ebnm)

set.seed(8) # For reproducibility of sampler results.

x <- c(28, 8, -3, 7, -1, 1, 18, 12) # Observations.
s <- c(15, 10, 16, 11, 9, 11, 10, 18) # Standard errors.
n <- length(x)

# By default, a posterior sampler isn't returned, so we set
#   output = output_all().
ebnm_res_mode0 <- ebnm::ebnm_point_normal(
  x = x, s = s, mode = 0, output = output_all()
)
mode0_samp <- ebnm_res_mode0$posterior_sampler(10000)
mode0_CI <- apply(mode0_samp, 2, quantile, c(.025, .975))

ebnm_res_estmode <- ebnm::ebnm_point_normal(
  x = x, s = s, mode = "estimate", output = output_all()
)
estmode_samp <- ebnm_res_estmode$posterior_sampler(10000)
estmode_CI <- apply(estmode_samp, 2, quantile, c(.025, .975))

modes <- paste("mode =", c("0", "\"estimate\""))
tib <- tibble(
  idx = factor(rep(LETTERS[1:n], 2)),
  x = rep(x, 2),
  s = rep(s, 2),
  theta_hat = c(
    ebnm_res_mode0$posterior$mean, ebnm_res_estmode$posterior$mean
  ),
  CI_lower = c(mode0_CI[1, ], estmode_CI[1, ]),
  CI_upper = c(mode0_CI[2, ], estmode_CI[2, ]),
  mode = factor(rep(modes, each = n), levels = modes)
)

ggplot(tib, aes(x = idx, y = x)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = x - 2 * s, ymax = x + 2 * s)) +
  geom_point(aes(x = idx, y = theta_hat), color = "red", size = 3) +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), color = "red") +
  theme_minimal() +
  labs(x = "\nSchool", y = "Treatment effect\n") +
  facet_wrap(~mode)

ggsave("../../figs/eight_schools.png", height = 3.5, width = 5)
