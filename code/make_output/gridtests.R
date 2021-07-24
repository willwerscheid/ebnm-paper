set.seed(666)
ntrials <- 20

if (exists("test") && test) {
  ns <- 10^(2:3)
} else {
  ns <- 10^(2:4)
}

maxcomp <- 5
families <- c("npmle", "normalmix", "unimix")

all_res <- tibble()
for (family in families) {
  for (n in ns) {
    if (family == "npmle") {
      n <- 3 * n / 10
    }
    cat("n =", n, "\n")
    for (i in 1:ntrials) {
      cat("  trial", i, "\n")
      ncomp <- sample(1:maxcomp)
      scale <- rexp(ncomp, rate = 1/5)
      prop <- runif(ncomp)
      prop <- prop / sum(prop)
      
      theta.scale <- sample(scale, n, replace = TRUE)
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
      for (K in c(30, 50, 80, 100, 120, 150, 200, 300)) {
        if (family == "npmle") {
          g_init <- ebnm:::init_g_for_npmle(x, s = 1, scale = sqrt(8 / n))
          ebnm_res <- do.call(ebnm_fn, list(x, s = 1, g_init = g_init, fix_g = FALSE))
        } else {
          g_scale <- ebnm:::default_scale(x, s = 1, mode = 0, min_K = K, max_K = K)
          ebnm_res <- do.call(ebnm_fn, list(x, s = 1, scale = g_scale))
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

saveRDS(all_res, "../../output/gridtests.rds")
