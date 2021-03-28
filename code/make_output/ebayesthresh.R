## Performance comparisons: ebnm_point_laplace vs. EBayesThresh.

cat("\nComparing ebnm_point_laplace and EBayesThresh.\n",
    "This takes about 20 minutes per test.\n",
    "Six tests are performed.\n\n")

null_sim <- function(n, s) {
  return(rnorm(n, sd = s))
}

pl_sim <- function(n, s) {
  pi0 <- rbeta(1, shape1 = 10, shape2 = 2)
  scale <- rgamma(1, shape = 4, rate = 1)
  
  theta <- rexp(n, rate = 1 / scale) * sample(c(-1, 1), n, replace = TRUE)
  theta <- theta * rbinom(n, 1, 1 - pi0)
  x <- theta + rnorm(n, sd = s)
  
  return(x)
}

pl_nzmode_sim <- function(n, s) {
  mode <- runif(1, -10, 10)
  x <- pl_sim(n, s) + mode
  
  return(x)
}

do_test <- function(n, sim_fn, homosked, nsim, mbtimes = 5L) {
  t <- rep(0, 2)
  ourbad <- 0
  theirbad <- 0
  theirfail <- 0

  for (i in 1:nsim) {
    if (homosked) {
      s <- 1
    } else {
      s <- sqrt(rexp(n))
    }
    
    x <- do.call(sim_fn, list(n = n, s = s))
    
    mb_res <- microbenchmark(
      ourres <- ebnm_point_laplace(x = x, s = s, optmethod = "nohess_nlm",
                                   control = list(gradtol = sqrt(.Machine$double.eps)),
                                   output = c("posterior_mean", "log_likelihood")),
      theirres <- tryCatch(
        ebayesthresh(x = x, sdev = s, a = NA, threshrule = "mean", 
                     universalthresh = FALSE, verbose = TRUE),
        error = function(x) list(w = NULL, a = NULL)
      ),
      times = mbtimes
    )
    
    ourllik <- ourres$log_likelihood
    if (is.null(theirres$w)) {
      theirllik <- -Inf
    } else {
      theirllik <- ebnm:::loglik_point_laplace(x, s, theirres$w, theirres$a, 0)
    }
    
    llik_diff <- ourllik - theirllik
    
    t <- t + summary(mb_res, unit = "s")$mean
    ourbad <- ourbad + 1L * (llik_diff < -n * sqrt(.Machine$double.eps))
    theirbad <- theirbad + 1L * (llik_diff > n * sqrt(.Machine$double.eps))
    theirfail <- theirfail + 1L * is.infinite(theirllik)
  }
  
  call <- match.call()
  res <- tibble(
    n = n,
    sim_fn = as.character(call$sim_fn),
    homosked = homosked,
    package = c("ebnm", "EBayesThresh"),
    mean_t = t / nsim,
    p_bad_soln = c(ourbad, theirbad) / nsim,
    p_fail = c(0, theirfail) / nsim
  )
  
  return(res)
}

set.seed(666)
all_res <- tibble()
if (exists("test") && test) {
  nsim <- 10^(2:3)
} else {
  ns <- 10^(3:5)
}

for (n in ns) {
  if (exists("test") && test) {
    nsim <- 1
  } else {
    nsim <- 10^6 / n
  }
  
  for (homosked in c(TRUE, FALSE)) {
    cat("  Testing n =", n, ", homosked =", homosked, "\n")
    all_res <- all_res %>%
      bind_rows(do_test(n = n, sim_fn = pl_sim, homosked = homosked, nsim = nsim)) %>%
      bind_rows(do_test(n = n, sim_fn = null_sim, homosked = homosked, nsim = nsim)) %>%
      bind_rows(do_test(n = n, sim_fn = pl_nzmode_sim, homosked = homosked, nsim = nsim)) 
  }
}

saveRDS(all_res, "../../output/ebayesthresh.rds")
