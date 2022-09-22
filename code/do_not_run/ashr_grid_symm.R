## Quality of ashr grid approximations for symmetric unimodal priors.
library(tidyverse)

log.add <- function(log.x, log.y) {
  c <- pmax(log.x, log.y)
  return(log(exp(log.x - c) + exp(log.y - c)) + c)
}

# assumes x > y:
log.minus <- function(log.x, log.y) {
  return(log(1 - exp(log.y - log.x)) + log.x)
}

UN.logdens <- function(x, a) {
  if (a == 0) {
    return(dnorm(x, log = TRUE))
  } else {
    return(log.minus(pnorm(a - x, log.p = TRUE), pnorm(-a - x, log.p = TRUE)) - log(2 * a))
  }
}

symmKLdiv <- function(a, omega, a.left, a.right) {
  KL.integrand <- function(x) {
    f.dens <- exp(UN.logdens(x, a = a))
    logf.dens <- UN.logdens(x, a = a)
    logg1.dens <- log(omega) + UN.logdens(x, a = a.left)
    logg2.dens <- log(1 - omega) + UN.logdens(x, a = a.right)
    logg.dens <- log.add(logg1.dens, logg2.dens)
    return(f.dens * (logf.dens - logg.dens))
  }
  return(integrate(KL.integrand, lower = -a - 5, upper = a + 5, rel.tol = sqrt(.Machine$double.eps))$value)
}

min.symmKLdiv <- function(a, a.left, a.right) {
  optres <- optimize(
    function(omega) symmKLdiv(a, omega, a.left, a.right), 
    interval = c(0, 1), 
    maximum = FALSE
  )
  return(optres)
}

ub.symmKLdiv <- function(a.left, a.right) {
  optres <- optimize(
    function(a) min.symmKLdiv(a, a.left, a.right)$objective, 
    interval = c(a.left, a.right), 
    maximum = TRUE
  )
  retlist <- list(
    KL.div = optres$objective,
    opt.a = optres$maximum,
    opt.omega = min.symmKLdiv(optres$maximum, a.left, a.right)$minimum
  )
  return(retlist)
}

find.next.gridpt <- function(a.left, targetKL, max.space, min.space) {
  uniroot.fn <- function(a.right) {
    return(ub.symmKLdiv(a.left, a.right)$KL.div - targetKL)
  }
  optres <- uniroot(uniroot.fn, c(a.left + min.space, a.left + max.space))
  return(optres$root)
}

print.progress <- function(K) {
  if (K %% 50 == 0) {
    cat("X\n")
  } else if (K %% 10 == 0) {
    cat("X")
  } else {
    cat(".")
  }
}

build.grid <- function(targetKL, max.K = 30, max.x = Inf, init.srch = c(0.1, 5), init.K = 5) {
  cat("Building grid ( KL =", targetKL, ")..")
  grid <- 0
  K <- 1
  
  is.space.increasing <- FALSE
  while (K < init.K && max(grid) < max.x) {
    next.gridpt <- find.next.gridpt(
      grid[K], targetKL, min.space = init.srch[1], max.space = init.srch[2]
    )
    grid <- c(grid, next.gridpt)
    K <- K + 1
    print.progress(K)
  }
  
  min.incr <- 1 + exp(1) * targetKL
  while (length(grid) < max.K && max(grid) < max.x) {
    last.space <- grid[length(grid)] - grid[length(grid) - 1]
    last.incr <- grid[length(grid)] / grid[length(grid) - 1]
    next.gridpt <- find.next.gridpt(
      grid[K], targetKL, min.space = last.space * min.incr, max.space = last.space * last.incr
    )
    grid <- c(grid, next.gridpt)
    K <- K + 1
    print.progress(K)
  }
  cat("\n")
  return(grid)
}

max.K <- 300
tib <- tibble()
for (log.KL in seq(-2, -2.75, by = -0.25)) {
  grid <- build.grid(10^(log.KL), max.K = max.K, init.K = 5, init.srch = c(1, 10))
  tib <- tib %>%
    bind_rows(tibble(KL = 10^(log.KL), idx = 1:length(grid), loc = grid))
}
for (log.KL in seq(-3, -4.75, by = -0.25)) {
  grid <- build.grid(10^(log.KL), max.K = max.K, init.K = 10, init.srch = c(0.2, 2))
  tib <- tib %>%
    bind_rows(tibble(KL = 10^(log.KL), idx = 1:length(grid), loc = grid))
}
for (log.KL in seq(-5, -5.75, by = -0.25)) {
  grid <- build.grid(10^(log.KL), max.K = max.K, init.K = 15, init.srch = c(0.1, 1))
  tib <- tib %>%
    bind_rows(tibble(KL = 10^(log.KL), idx = 1:length(grid), loc = grid))
}
for (log.KL in seq(-6, -6.75, by = -0.25)) {
  grid <- build.grid(10^(log.KL), max.K = max.K, init.K = 25, init.srch = c(0.05, 1))
  tib <- tib %>%
    bind_rows(tibble(KL = 10^(log.KL), idx = 1:length(grid), loc = grid))
}
for (log.KL in seq(-7, -8, by = -0.25)) {
  grid <- build.grid(10^(log.KL), max.K = max.K, init.K = 50, init.srch = c(0.01, 0.5))
  tib <- tib %>%
    bind_rows(tibble(KL = 10^(log.KL), idx = 1:length(grid), loc = grid))
}
