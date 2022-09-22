## Estimating difference in objective between ebnm solutions and 'optimal' solutions.

set.seed(666)
ntrials <- 20

if (exists("test") && test) {
  ns <- 10^(1:2)
  Ks <- c(30, 50)
} else {
  ns <- c(30, 100, 300, 1000)
  Ks <- c(30, 50, 80, 100, 120, 150, 200, 300)
}

cat(paste0("\nEstimating quality of grids for different prior families.\n",
           "  n ranges from ", min(ns), " to ", max(ns), ".\n\n"))

maxcomp <- 5
families <- c("normalmix", "unimix", "npmle")

t_begin <- Sys.time()
all_res <- tibble()
for (family in families) {
  cat("  Prior family:", family, "\n")
  for (n in ns) {
    cat("    n:", n, "\n")
    for (i in 1:ntrials) {
      cat("      trial", i, "\n")
      ncomp <- sample(maxcomp, 1)
      scale <- rexp(ncomp, rate = 1/5)
      
      if (ncomp == 1) {
        theta.scale <- rep(scale, n)
      } else {
        prop <- runif(ncomp)
        prop <- prop / sum(prop)
        theta.scale <- sample(scale, n, replace = TRUE, prob = prop)
      }
      
      if (family == "normalmix") {
        theta <- rnorm(n, sd = theta.scale)
        ebnm_fn <- ebnm_normal_scale_mixture
      } else if (family == "unimix") {
        theta <- runif(n, min = -theta.scale, max = theta.scale)
        ebnm_fn <- ebnm_unimodal_symmetric
      } else if (family == "npmle") {
        theta <- theta.scale
        ebnm_fn <- ebnm_npmle
      }
      x <- theta + rnorm(n)

      ebnm_def <- do.call(ebnm_fn, list(x = x, s = 1))

      res <- tibble()
      for (K in Ks) {
        if (family == "npmle") {
          g_init <- ebnm:::init_g_for_npmle(x, s = 1, scale = (max(x) - min(x)) / (K - 1))
          ebnm_res <- do.call(ebnm_fn, list(x, s = 1, g_init = g_init, fix_g = FALSE))
        } else if (family == "normalmix") {
          g_scale <- ebnm:::default_smn_scale(x, s = 1, mode = 0, min_K = K, max_K = K)
          ebnm_res <- do.call(ebnm_fn, list(x, s = 1, scale = g_scale))
        } else if (family == "unimix") {
          g_scale <- ebnm:::default_symmuni_scale(x, s = 1, mode = 0, min_K = K, max_K = K)
          ebnm_res <- do.call(ebnm_fn, list(x, s = 1, scale = g_scale))
        } else {
          stop("Unrecognized prior family.")
        }

        res <- res %>%
          bind_rows(tibble(K = K, llik = ebnm_res$log_likelihood))
      }

      res <- res %>%
        mutate(llik_diff = llik - ebnm_def$log_likelihood)

      all_res <- all_res %>%
        bind_rows(tibble(
          family = family,
          trial = i,
          n = n,
          diff = res %>% summarize(max_diff = max(llik_diff)) %>% pull(max_diff)
        ))
    }
  }
}

t_elapsed <- Sys.time() - t_begin
cat("  Done. Time elasped:",
    round(as.numeric(t_elapsed, units = "mins"), 1), "minutes.\n")

saveRDS(all_res, "../../output/gridtests.rds")
rm(res)
