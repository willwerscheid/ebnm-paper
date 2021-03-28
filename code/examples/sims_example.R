## Simulations examples.

cat("\nSimulation examples...\n\n")

# 1. Mixture of point mass at zero and t_5
point_t <- function(n) {
  samp <- 1.5 * rt(n, df = 5) 
  samp <- samp * sample(c(0, 1), n, replace = TRUE, prob = c(0.5, 0.5))
  return(samp + rnorm(n))
}

# Simulate data and estimate posterior means
set.seed(1)
n <- 1000
x <- point_t(n)
ebnm_fns <- c("ebnm_normal_scale_mixture",
              "ebnm_point_normal",
              "ebnm_horseshoe")
ebnm_res <- list()
for (ebnm_fn in ebnm_fns) {
  ebnm_res[[ebnm_fn]] <- do.call(ebnm_fn, list(x = x, s = 1))
}

# Compare log likelihoods
cat("Point-t log likelihoods:\n")
print(sapply(ebnm_res, `[[`, "log_likelihood"))
cat("\n")

# Plot posterior means
tib <- tibble()
for (ebnm_fn in ebnm_fns) {
  tib <- tib %>%
    bind_rows(tibble(fn = ebnm_fn,
                     x = x,
                     postmean = ebnm_res[[ebnm_fn]]$posterior$mean))
}
ggplot(tib, aes(x = x, y = postmean, col = fn)) + 
  geom_point(size = 0.5) +
  geom_abline(slope = 1, linetype = "dashed") +
  labs(x = "x", y = "posterior mean", col = "ebnm function") +
  theme_minimal()

ggsave("../../figs/pointt.png", width = 6, height = 3.375)


# 2. Mixture of point mass at zero and Unif[-5, 10]
asymm_tophat <- function(n) {
  samp <- runif(n, -5, 10) 
  samp <- samp * sample(c(0, 1), n, replace = TRUE, prob = c(0.5, 0.5))
  return(samp + rnorm(n))
}

set.seed(1)
x <- asymm_tophat(n)
ebnm_fns <- c("ebnm_unimodal",
              "ebnm_unimodal_symmetric",
              "ebnm_normal_scale_mixture")
ebnm_res <- list()
for (ebnm_fn in ebnm_fns) {
  ebnm_res[[ebnm_fn]] <- do.call(ebnm_fn, list(x = x, s = 1))
}

cat("Asymmetric tophat log likelihoods:\n")
print(sapply(ebnm_res, `[[`, "log_likelihood"))
cat("\n")

# Plot posterior means
tib <- tibble()
for (ebnm_fn in ebnm_fns) {
  tib <- tib %>%
    bind_rows(tibble(fn = ebnm_fn,
                     x = x,
                     postmean = ebnm_res[[ebnm_fn]]$posterior$mean))
}
ggplot(tib, aes(x = x, y = postmean, col = fn)) + 
  geom_point(size = 0.5) +
  geom_abline(slope = 1, linetype = "dashed") +
  labs(x = "x", y = "posterior mean", col = "ebnm function") +
  theme_minimal()

ggsave("../../figs/tophat.png", width = 6, height = 3.375)


# 3. Mixture of point mass at zero and point mass at 5
two_pm <- function(n) {
  samp <- sample(c(0, 5), n, replace = TRUE, prob = c(0.8, 0.2))
  return(samp + rnorm(n))
}

set.seed(1)
x <- two_pm(n)
ebnm_fns <- c("ebnm_npmle",
              "ebnm_deconvolver",
              "ebnm_unimodal")
ebnm_res <- list()
for (ebnm_fn in ebnm_fns) {
  ebnm_res[[ebnm_fn]] <- do.call(ebnm_fn, list(x = x, s = 1))
}

cat("Two-pointmass log likelihoods:\n")
print(sapply(ebnm_res, `[[`, "log_likelihood"))
cat("\n")

# Plot posterior means
tib <- tibble()
for (ebnm_fn in ebnm_fns) {
  tib <- tib %>%
    bind_rows(tibble(fn = ebnm_fn,
                     x = x,
                     postmean = ebnm_res[[ebnm_fn]]$posterior$mean))
}
ggplot(tib, aes(x = x, y = postmean, col = fn)) + 
  geom_point(size = 0.5) +
  geom_abline(slope = 1, linetype = "dashed") +
  labs(x = "x", y = "posterior mean", col = "ebnm function") +
  theme_minimal()

ggsave("../../figs/twopm.png", width = 6, height = 3.375)
