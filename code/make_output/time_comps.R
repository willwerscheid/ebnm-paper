## Timing different prior families.

cat("\nTiming different prior families.\n",
    "n ranges from 1e3 to 1e6.\n\n")

# Mixture of point mass at zero and t_5
point_t <- function(n) {
  samp <- 1.5 * rt(n, df = 5)
  samp <- samp * sample(c(0, 1), n, replace = TRUE, prob = c(0.5, 0.5))
  return(samp + rnorm(n))
}

if (exists("test") && test) {
  ns <- 10^seq(2, 3, by = 0.5)
  mb_times <- rep(2, length(ns))
} else {
  ns <- 10^seq(3, 6, by = 0.5)
  mb_times <- c(20, 20, 10, 10, 3, 1, 1)
}

ebnm_fns <- c("ebnm_normal",
              "ebnm_point_normal",
              "ebnm_point_laplace",
              "ebnm_normal_scale_mixture",
              "ebnm_unimodal_symmetric",
              "ebnm_unimodal",
              "ebnm_deconvolver",
              "ebnm_npmle",
              "ebnm_horseshoe")

set.seed(1)
res <- tibble()
for (i in 1:length(ns)) {
  n <- ns[i]
  cat("  n:", n, "\n")
  
  x <- point_t(n)
  
  test_fns <- ebnm_fns
  if (n == 10^6) {
    test_fns <- setdiff(ebnm_fns, c("ebnm_npmle", "ebnm_horseshoe"))
  } else {
    test_fns <- ebnm_fns
  }
  
  mb_tests <- lapply(test_fns, function(fn) {
    bquote(do.call(.(fn), list(x = x, s = 1)))
  })
  
  mb_res <- microbenchmark::microbenchmark(
    list = mb_tests,
    times = mb_times[i]
  )
  
  res <- res %>%
    bind_rows(tibble(
      n = n, 
      fn = test_fns, 
      t = summary(mb_res, unit = "s")$mean
    ))
}

cat("  Done.\n")

res <- res %>%
  mutate(fn = fct_relevel(fn, rev(ebnm_fns)))

saveRDS(res, "../../output/timecomps.rds")
rm(ebnm_res)
