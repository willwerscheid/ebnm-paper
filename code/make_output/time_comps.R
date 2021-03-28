## Timing different prior families.

cat("\nTiming different prior families.\n",
    "This takes about an hour.\n",
    "n ranges from 1e3 to 1e6.\n\n")

# Mixture of point mass at zero and t_5
point_t <- function(n) {
  samp <- 1.5 * rt(n, df = 5) 
  samp <- samp * sample(c(0, 1), n, replace = TRUE, prob = c(0.5, 0.5))
  return(samp + rnorm(n))
}

if (exists("test") && test) {
  ns <- 10^seq(2, 3, by = 0.5)
} else {
  ns <- 10^seq(3, 6, by = 0.5)
}

ebnm_fns <- c("ebnm_point_normal",
              "ebnm_point_laplace",
              "ebnm_normal_scale_mixture",
              "ebnm_unimodal_symmetric",
              "ebnm_unimodal",
              "ebnm_deconvolver",
              "ebnm_npmle",
              "ebnm_horseshoe")

set.seed(1)
res <- tibble()
for (n in ns) {
  cat("  n:", n, "\n")
  x <- point_t(n)
  for (ebnm_fn in ebnm_fns) {
    cat("    ebnm fn:", ebnm_fn, "\n")
    t <- system.time(ebnm_res <- do.call(ebnm_fn, list(x = x, s = 1)))[3]
    res <- res %>%
      bind_rows(tibble(n = n, fn = ebnm_fn, t = t, llik = ebnm_res$log_likelihood))
  }
}

res <- res %>%
  mutate(fn = fct_relevel(fn, rev(ebnm_fns)))

saveRDS(res, "../../output/timecomps.rds")
