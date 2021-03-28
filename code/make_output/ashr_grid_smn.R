## Quality of ashr grid approximations for scale mixtures of normals.

cat("\nEstimating quality of ashr grid approximations for scale mixtures of",
    "normals.\n This takes roughly 30 seconds per value of m.\n",
    "Values of m range from 1.1 to 2.\n\n")

set.seed(666)

# To sample from different distributions N(0, s2), use a single sample from 
#   N(0, 1) and scale. This method ensures a smooth optimization objective.
smnKLdiv <- function(s2, omega, m, samp) {
  samp <- samp * sqrt(s2)
  
  tmp <- omega * exp(0.5 * samp^2 * (1 / s2 - 1)) 
  tmp <- tmp + (1 - omega) * exp(0.5 * samp^2 * (1 / s2 - 1 / m)) / sqrt(m)
  tmp <- -log(tmp) - 0.5 * log(s2)
  
  return(mean(tmp))
}

min_smnKLdiv <- function(s2, m, samp) {
  optres <- optimize(
    function(omega) smnKLdiv(s2, omega, m, samp), 
    interval = c(0, 1), 
    maximum = FALSE
  )
  return(optres$objective)
}

ub_smnKLdiv <- function(m, samp) {
  optres <- optimize(
    function(s2) min_smnKLdiv(s2, m, samp), 
    interval = c(1, m), 
    maximum = TRUE
  )
  return(optres$objective)
}

if (exists("test") && test) {
  m <- seq(1.1, 2, by = 0.1)
} else {
  m <- seq(1.1, 2, by = 0.025)
}

if (exists("test") && test) {
  sampsize <- 100000
} else {
  sampsize <- 10000000
}
samp <- rnorm(sampsize)

ub <- numeric(length(m))
for (i in 1:length(m)) {
  cat("  m:", m[i], "\n")
  ub[i] <- ub_smnKLdiv(m[i], samp)
}

tib <- tibble(m = m, ub = ub)

saveRDS(tib, "../../output/smngrid.rds")
