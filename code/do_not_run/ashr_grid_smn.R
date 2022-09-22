## Quality of ashr grid approximations for scale mixtures of normals.

log.add <- function(log.x, log.y) {
  c <- pmax(log.x, log.y)
  return(log(exp(log.x - c) + exp(log.y - c)) + c)
}

smnKLdiv <- function(s2, omega, m) {
  KL.integrand <- function(x) {
    f.dens <- dnorm(x, mean = 0, sd = sqrt(s2))
    logf.dens <- dnorm(x, mean = 0, sd = sqrt(s2), log = TRUE)
    logg1.dens <- log(omega) + dnorm(x, mean = 0, sd = 1, log = TRUE)
    logg2.dens <- log(1 - omega) + dnorm(x, mean = 0, sd = sqrt(m), log = TRUE)
    logg.dens <- log.add(logg1.dens, logg2.dens)
    return(f.dens * (logf.dens - logg.dens))
  }
  return(integrate(KL.integrand, -Inf, Inf)$value)
}

min.smnKLdiv <- function(s2, m) {
  optres <- optimize(
    function(omega) smnKLdiv(s2, omega, m), 
    interval = c(0, 1), 
    maximum = FALSE
  )
  return(optres)
}

ub.smnKLdiv <- function(m) {
  optres <- optimize(
    function(s2) min.smnKLdiv(s2, m)$objective, 
    interval = c(1, m), 
    maximum = TRUE
  )
  retlist <- list(
    KL.div = optres$objective,
    opt.s2 = optres$maximum,
    opt.omega = min.smnKLdiv(optres$maximum, m)$minimum
  )
  return(retlist)
}

m <- exp(exp(seq(log(log(1.025)), log(log(4)), length.out = 100)))
ub <- numeric(length(m))

for (i in 1:length(m)) {
  cat("  m:", m[i], "\n")
  ub[i] <- ub.smnKLdiv(m[i])$KL.div
}
