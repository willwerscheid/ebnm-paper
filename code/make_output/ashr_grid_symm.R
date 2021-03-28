## Quality of ashr grid approximations for symmetric unimodal priors.

cat("\nEstimating quality of ashr grid approximations for symmetric unimodal",
    "priors.\n This takes roughly five minutes per function evaluation.\n",
    "KL divergences range from 1e-2 to 1e-6.\n\n")

set.seed(666)

lik_UN <- function(a, samp) {
  return((pnorm(a - samp) - pnorm(-a - samp)) / (2 * a))
}

# To sample from different distributions Uniform[-a, a] + N(0, 1), use a single sample 
#   from Unif[-1, 1] and scale. This method ensures a smooth optimization objective.
symmKLdiv <- function(a, omega, aleft, aright, unif_samp, norm_samp) {
  samp <- a * unif_samp + norm_samp
  
  tmp <- log(lik_UN(a, samp))
  tmp <- tmp - log(omega * lik_UN(aleft, samp) + (1 - omega) * lik_UN(aright, samp))
  
  return(mean(tmp))
}

min_symmKLdiv <- function(a, aleft, aright, unif_samp, norm_samp) {
  optres <- optimize(
    function(omega) symmKLdiv(a, omega, aleft, aright, unif_samp, norm_samp),
    interval = c(0, 1), 
    maximum = FALSE
  )
  return(optres$objective)
}

ub_symmKLdiv <- function(aleft, aright, unif_samp, norm_samp) {
  optres <- optimize(
    function(a) min_symmKLdiv(a, aleft, aright, unif_samp, norm_samp),
    interval = c(aleft, aright), 
    maximum = TRUE
  )
  return(optres$objective)
}

find_next_gridpt <- function(aleft, targetKL, unif_samp, norm_samp, max_space) {
  uniroot_fn <- function(aright) {
    ub_symmKLdiv(aleft, aright, unif_samp, norm_samp) - targetKL
  }
  optres <- uniroot(uniroot_fn, c(aleft + 1e-6, aleft + max_space))
  return(optres$root)
}

find_min_space <- function(targetKL, unif_samp, norm_samp, max_space) {
  optres <- optimize(
    function(aleft) {
      cat("    Evaluating at", aleft, "\n")
      aright <- find_next_gridpt(aleft, targetKL, unif_samp, norm_samp, max_space)
      cat("      Result:", aright - aleft, "\n")
      return(aright - aleft)
    },
    interval = c(1e-6, 5),
    tol = 0.01,
    maximum = FALSE
  )
  return(optres$objective)
}

# The sample size needs to be relatively small.
if (exists("test") && test) {
  sampsize <- 1000
} else {
  sampsize <- 1000000
}
unif_samp <- runif(sampsize, min = -1, max = 1)
norm_samp <- rnorm(sampsize)

tib <- NULL
if (exists("test") && test) {
  KLs <- 10^(seq(-2, -3, by = -0.5))
} else {
  KLs <- 10^(seq(-2, -6, by = -0.25))
}

for (KL in KLs) {
  cat("  KL:", KL, "\n")
  min_space <- find_min_space(KL, unif_samp, norm_samp, max_space = 5)
  if (is.null(tib)) {
    tib <- tibble(KL = KL, min_space = min_space)
  } else {
    tib <- tib %>%
      add_case(KL = KL, min_space = min_space)
  }
}

saveRDS(tib, "../../output/symmgrid.rds")
