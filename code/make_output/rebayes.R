## Performance comparisons: ebnm_npmle vs. REBayes.

cat("\nComparing ebnm_npmle and REBayes.\n",
    "24 tests are performed.\n\n")

pl_sim <- function(n, s) {
  pi0 <- rbeta(1, shape1 = 10, shape2 = 2)
  scale <- rgamma(1, shape = 4, rate = 1)
  
  theta <- rexp(n, rate = 1 / scale) * sample(c(-1, 1), n, replace = TRUE)
  theta <- theta * rbinom(n, 1, 1 - pi0)
  x <- theta + rnorm(n, sd = s)
  
  return(x)
}

do_test <- function(n, homosked, nsim, n_gridpts, mbtimes = 5L) {
  t <- rep(0, 2)
  llik_diff <- 0
  theirfail <- 0

  for (i in 1:nsim) {
    if (homosked) {
      s <- 1
    } else {
      s <- sqrt(rexp(n))
    }
    
    x <- pl_sim(n = n, s = s)
    
    scale <- (max(x) - min(x)) / (n_gridpts - 1)
    
    mb_res <- microbenchmark(
      ourres <- ebnm_npmle(x = x, s = s, scale = scale,
                           control = list(convtol.sqp = 1e-4),
                           output = c("posterior_mean", "log_likelihood")),
      theirres <- GLmix(x, sigma = s, v = n_gridpts, control = list(rtol = 1e-4)),
      times = mbtimes
    )
    
    ourllik <- ourres$log_likelihood
    theirllik <- theirres$logLik
    
    if (theirres$status == "OPTIMAL") {
      llik_diff <- llik_diff + ourllik - theirllik
      t <- t + summary(mb_res, unit = "s")$mean
    } else {
      theirfail <- theirfail + 1
    }
  }
  
  call <- match.call()
  res <- tibble(
    n = n,
    n_gridpts = n_gridpts,
    homosked = homosked,
    package = c("ebnm", "REBayes"),
    mean_t = t / nsim,
    mean_llik_diff = c(1, -1) * llik_diff / nsim,
    prop_fail = c(0, theirfail / nsim)
  )
  
  return(res)
}

set.seed(666)
all_res <- tibble()
if (exists("test") && test) {
  ns <- 10^(2:3)
  ngridpts <- c(10, 30)
} else {
  ns <- 10^(3:5)
  ngridpts <- c(10, 30, 100, 300)
}

for (n in ns) {
  for (ngridpt in ngridpts) {
    if (exists("test") && test) {
      nsim <- 10
    } else {
      nsim <- min(100, ceiling(max(ngridpts) * max(ns) / (ngridpt * n)))
    }
    
    for (homosked in c(TRUE, FALSE)) {
      cat("  Testing n =", n, ", n_gridpts =", ngridpt, ", homosked =", homosked, "\n")
      all_res <- all_res %>%
        bind_rows(do_test(n = n, homosked = homosked, nsim = nsim, n_gridpts = ngridpt))
    }
  }
}

saveRDS(all_res, "../../output/rebayes.rds")
