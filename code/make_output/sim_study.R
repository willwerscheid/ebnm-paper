cat("\nSimulation study...\n\n")

point_normal <- function(n) {
  samp <- rnorm(n, sd = 2) 
  samp <- samp * sample(c(0, 1), n, replace = TRUE, prob = c(0.9, 0.1))
  return(list(theta = samp, x = samp + rnorm(n)))
}

point_t <- function(n) {
  samp <- 1.5 * rt(n, df = 5) 
  samp <- samp * sample(c(0, 1), n, replace = TRUE, prob = c(0.8, 0.2))
  return(list(theta = samp, x = samp + rnorm(n)))
}

asymm_tophat <- function(n) {
  samp <- runif(n, -5, 10) 
  samp <- samp * sample(c(0, 1), n, replace = TRUE, prob = c(0.5, 0.5))
  return(list(theta = samp, x = samp + rnorm(n)))
}

if (exists("test") && test) {
  nsim <- 1
  n <- 100
} else {
  nsim <- 10
  n <- 1000
}

ebnm_fns <- c("ebnm_normal",
              "ebnm_point_normal",
              "ebnm_point_laplace",
              "ebnm_normal_scale_mixture",
              "ebnm_unimodal_symmetric",
              "ebnm_unimodal",
              "ebnm_npmle",
              "ebnm_horseshoe",
              "ebnm_deconvolver")
output <- c("posterior_mean",
            "log_likelihood",
            "posterior_sampler")
sim_fns <- c("point_normal",
             "point_t",
             "asymm_tophat")
             
all_res <- tibble()
for (sim_fn in sim_fns) {
  cat("  Sim function:", sim_fn, "\n")
  for (i in 1:nsim) {
    cat("    trial", i, "\n")
    set.seed(i)
    sim_data <- do.call(sim_fn, list(n = n))
    ebnm_res <- list()
    for (ebnm_fn in ebnm_fns) {
      ebnm_res[[ebnm_fn]] <- do.call(ebnm_fn, list(x = sim_data$x, s = 1, output = output))
    }
    
    llik <- sapply(ebnm_res, `[[`, "log_likelihood")
    llik <- llik - max(llik)
    
    rmse <- sapply(ebnm_res, function(res) {
      return(sqrt(mean((res$posterior$mean - sim_data$theta)^2)))
    })
    
    confint_cov <- sapply(ebnm_res, function(res) {
      zz <- capture.output({
        samp <- res$posterior_sampler(10000)
      })
      ci <- apply(samp, 2, quantile, probs = c(0.05, 0.95))
      return(1 - mean(sim_data$theta < ci[1, ] | sim_data$theta > ci[2, ]))
    })
    
    tib <- tibble(
      SimFn = sim_fn,
      Function = names(llik), 
      LogLikelihood = llik, 
      RMSE = rmse, 
      ConfIntCov = confint_cov,
      SimNumber = i
    )
    
    all_res <- all_res %>%
      bind_rows(tib)
  }
}

saveRDS(all_res, "../../output/simstudy.rds")
rm(all_res)
