## Performance comparisons: ebnm_normal_scale_mixture vs. ashr.

cat("\nComparing ebnm_normal_scale_mixture and ashr.\n",
    "This takes about 30 minutes per test.\n",
    "Six tests are performed.\n\n")

null_sim <- function(n, s) {
  return(rnorm(n, sd = s))
}

pt_sim <- function(n, s) {
  pi0 <- rbeta(1, shape1 = 10, shape2 = 2)
  scale <- rgamma(1, shape = 4, rate = 1)
  
  theta <- scale * rt(n, df = 5)
  theta <- theta * rbinom(n, 1, 1 - pi0)
  x <- theta + rnorm(n, sd = s)
  
  return(x)
}

pt_nzmode_sim <- function(n, s) {
  mode <- runif(1, -10, 10)
  x <- pt_sim(n, s) + mode
  
  return(x)
}

do_test <- function(n, sim_fn, homosked, nsim, mbtimes = 5L) {
  t <- rep(0, 2)

  for (i in 1:nsim) {
    if (homosked) {
      s <- 1
    } else {
      s <- sqrt(rexp(n))
    }
    
    x <- do.call(sim_fn, list(n = n, s = s))
    
    mb_res <- microbenchmark(
      ourres <- ebnm_normal_scale_mixture(
        x = x, 
        s = s, 
        output = c("fitted_g", "posterior_mean", "posterior_sd", "log_likelihood")
      ),
      theirres <- ash(
        betahat = x, 
        sebetahat = s, 
        mixcompdist = "normal", 
        prior = "uniform", 
        outputlevel = 2
      ),
      times = mbtimes
    )
    
    t <- t + summary(mb_res, unit = "s")$mean
  }
  
  call <- match.call()
  res <- tibble(
    n = n,
    sim_fn = as.character(call$sim_fn),
    homosked = homosked,
    package = c("ebnm", "ashr"),
    mean_t = t / nsim,
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
      bind_rows(do_test(n = n, sim_fn = pt_sim, homosked = homosked, nsim = nsim)) %>%
      bind_rows(do_test(n = n, sim_fn = null_sim, homosked = homosked, nsim = nsim)) %>%
      bind_rows(do_test(n = n, sim_fn = pt_nzmode_sim, homosked = homosked, nsim = nsim)) 
  }
}

saveRDS(all_res, "../../output/ebnmvashr.rds")
