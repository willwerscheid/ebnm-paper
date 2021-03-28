## Optimal ashr grids for symmetric unimodal priors.

cat("\nConstructing 'optimal' ashr grids for symmetric unimodal priors.\n",
    "This takes roughly five minutes per gridpoint.\n",
    "Grids are constructed on [0, 10].\n\n")

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

build_grid <- function(targetKL, unif_samp, norm_samp, max_space, max_val) {
  next_gridpt <- 1e-6
  grid <- next_gridpt
  while (next_gridpt < max_val) {
    next_gridpt <- find_next_gridpt(next_gridpt, targetKL, unif_samp, norm_samp, max_space)
    cat("    Gridpoint at:", next_gridpt, "\n")
    grid <- c(grid, next_gridpt)
  }
  return(grid[-length(grid)])
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
  KLs <- 10^(seq(-2, -3, by = -1))
} else {
  KLs <- 10^(seq(-2, -6, by = -1))
}

for (KL in KLs) {
  cat("  KL:", KL, "\n")
  grid <- build_grid(KL, unif_samp, norm_samp, max_space = 5, max_val = 10)
  if (is.null(tib)) {
    tib <- tibble(KL = KL, idx = 1:length(grid), grid = grid)
  } else {
    tib <- tib %>%
      bind_rows(tibble(KL = KL, idx = 1:length(grid), grid = grid))
  }
}

saveRDS(tib, "../../output/optgrids.rds")
