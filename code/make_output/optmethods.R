## Performance comparisons: ebnm optmethods.

cat("\nComparing ebnm optimization methods.\n",
    "This takes about 20 minutes per test.\n",
    "24 tests are performed.\n\n")

null_sim <- function(n, s) {
  return(rnorm(n, sd = s))
}

pn_sim <- function(n, s) {
  pi0 <- rbeta(1, shape1 = 10, shape2 = 2)
  scale <- rgamma(1, shape = 4, rate = 1)
  
  theta <- rnorm(n, sd = scale)
  theta <- theta * rbinom(n, 1, 1 - pi0)
  x <- theta + rnorm(n, sd = s)
  
  return(x)
}

pn_nzmode_sim <- function(n, s) {
  mode <- runif(1, -10, 10)
  x <- pn_sim(n, s) + mode
  
  return(x)
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

ebnm_pn_estmode <- function(...) {
  return(ebnm_point_normal(..., mode = "estimate"))
}

ebnm_pl_estmode <- function(...) {
  return(ebnm_point_laplace(..., mode =))
}

do_test <- function(n, sim_fn, homosked, ebnm_fn, nsim, mbtimes = 5L) {
  t <- rep(0, 6)
  bad_soln <- rep(0, 6)
  fail <- rep(0, 6)
  
  for (i in 1:nsim) {
    if (homosked) {
      s <- 1
    } else {
      s <- sqrt(rexp(n))
    }
    
    x <- do.call(sim_fn, list(n = n, s = s))
    
    ebnm_args <- list(x = x, s = s, output = "log_likelihood")
    
    llik <- numeric(6)
    mb_res <- microbenchmark(
      llik[1] <- unlist(
        do.call(ebnm_fn, 
                c(ebnm_args, 
                  list(optmethod = "nlm", 
                       control = list(gradtol = sqrt(.Machine$double.eps)))))
      ),
      llik[2] <- unlist(
        do.call(ebnm_fn, 
                c(ebnm_args, 
                  list(optmethod = "nohess_nlm", 
                       control = list(gradtol = sqrt(.Machine$double.eps)))))
      ),
      llik[3] <- unlist(
        do.call(ebnm_fn, 
                c(ebnm_args, 
                  list(optmethod = "nograd_nlm", 
                       control = list(gradtol = sqrt(.Machine$double.eps)))))
      ),
      llik[4] <- tryCatch(
        unlist(
          do.call(ebnm_fn,
                  c(ebnm_args,
                    list(optmethod = "lbfgsb",
                         control = list(pgtol = sqrt(.Machine$double.eps)))))
        ),
        error = function(x) -Inf
      ),
      llik[5] <- tryCatch(
        unlist(
          do.call(ebnm_fn,
                  c(ebnm_args,
                    list(optmethod = "nograd_lbfgsb",
                         control = list(pgtol = sqrt(.Machine$double.eps)))))
        ),
        error = function(x) -Inf
      ),
      llik[6] <- unlist(
        do.call(ebnm_fn, 
                c(ebnm_args,
                  list(optmethod = "trust",
                       control = list(fterm = sqrt(.Machine$double.eps)))))
      ),
      times = mbtimes
    )
    
    llik <- llik - max(llik)
    t <- t + summary(mb_res, unit = "s")$mean
    bad_soln <- bad_soln + 1L * (llik < -n * sqrt(.Machine$double.eps))
    fail <- fail + 1L * is.infinite(llik)
  }
  
  call <- match.call()
  res <- tibble(
    n = n,
    ebnm_fn = as.character(call$ebnm_fn),
    sim_fn = as.character(call$sim_fn),
    homosked = homosked,
    optmethod = c("nlm", "nohess_nlm", "nograd_nlm", "lbfgsb", "nograd_lbfgsb", "trust"),
    mean_t = t / nsim,
    p_bad_soln = bad_soln / nsim,
    p_fail = fail / nsim
  )
  
  return(res)
}

set.seed(666)
all_res <- tibble()
if (exists("test") && test) {
  ns <- 10^(2:3)
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
    cat("  Testing ebnm_point_normal: n =", n, ", homosked =", homosked, "\n")
    all_res <- all_res %>%
      bind_rows(do_test(n = n, sim_fn = pn_sim, homosked = homosked, 
                        ebnm_fn = ebnm_point_normal, nsim = nsim)) %>%
      bind_rows(do_test(n = n, sim_fn = null_sim, homosked = homosked, 
                        ebnm_fn = ebnm_point_normal, nsim = nsim)) %>%
      bind_rows(do_test(n = n, sim_fn = pn_nzmode_sim, homosked = homosked, 
                        ebnm_fn = ebnm_point_normal, nsim = nsim)) 
    cat("  Testing ebnm_pn_estmode: n =", n, ", homosked = ", homosked, "\n")
    all_res <- all_res %>%
      bind_rows(do_test(n = n, sim_fn = pn_nzmode_sim, homosked = homosked, 
                        ebnm_fn = ebnm_pn_estmode, nsim = nsim)) %>%
      bind_rows(do_test(n = n, sim_fn = null_sim, homosked = homosked, 
                        ebnm_fn = ebnm_pn_estmode, nsim = nsim)) %>%
      bind_rows(do_test(n = n, sim_fn = pt_nzmode_sim, homosked = homosked, 
                        ebnm_fn = ebnm_pn_estmode, nsim = nsim))
    cat("  Testing ebnm_point_laplace: n =", n, ", homosked =", homosked, "\n")
    all_res <- all_res %>%
      bind_rows(do_test(n = n, sim_fn = pl_sim, homosked = homosked, 
                        ebnm_fn = ebnm_point_laplace, nsim = nsim)) %>%
      bind_rows(do_test(n = n, sim_fn = null_sim, homosked = homosked, 
                        ebnm_fn = ebnm_point_laplace, nsim = nsim)) %>%
      bind_rows(do_test(n = n, sim_fn = pl_nzmode_sim, homosked = homosked, 
                        ebnm_fn = ebnm_point_laplace, nsim = nsim)) 
    cat("  Testing ebnm_pl_estmode: n =", n, ", homosked = ", homosked, "\n")
    all_res <- all_res %>%
      bind_rows(do_test(n = n, sim_fn = pl_nzmode_sim, homosked = homosked, 
                        ebnm_fn = ebnm_pl_estmode, nsim = nsim)) %>%
      bind_rows(do_test(n = n, sim_fn = null_sim, homosked = homosked, 
                        ebnm_fn = ebnm_pl_estmode, nsim = nsim)) %>%
      bind_rows(do_test(n = n, sim_fn = pt_nzmode_sim, homosked = homosked, 
                        ebnm_fn = ebnm_pl_estmode, nsim = nsim))
  }
}

saveRDS(all_res, "../../output/optmethods.rds")
