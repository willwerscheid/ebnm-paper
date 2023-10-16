###### SIMULATION STUDY -----------------------------------------------

cat("\nSIMULATION STUDY...\n\n")

library(tidyverse)
library(ebnm)
library(gt)

#### Simulation functions ----

normal <- function(n) {
  samp <- rnorm(n, sd = 2) 
  return(list(theta = samp, x = samp + rnorm(n)))
}

point_t <- function(n) {
  samp <- 1.5 * rt(n, df = 5) 
  samp <- samp * sample(c(0, 1), n, replace = TRUE, prob = c(0.8, 0.2))
  return(list(theta = samp, x = samp + rnorm(n)))
}

asymm_tophat <- function(n) {
  samp <- runif(n, -5, 10) 
  return(list(theta = samp, x = samp + rnorm(n)))
}

#### Run simulations ----

set.seed(666)

if (exists("test") && test) {
  nsim <- 1
  n <- 100
} else {
  nsim <- 10
  n <- 1000
}

ebnm_fns <- c("ebnm_flat",
              "ebnm_normal",
              "ebnm_point_normal",
              "ebnm_point_laplace",
              "ebnm_normal_scale_mixture",
              "ebnm_unimodal_symmetric",
              "ebnm_unimodal",
              "ebnm_npmle",
              "ebnm_deconvolver",
              "ebnm_horseshoe")
output <- c("posterior_mean",
            "log_likelihood",
            "posterior_sampler")
sim_fns <- c("normal",
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
    
    llik <- sapply(ebnm_res, logLik)
    llik <- llik - max(llik, na.rm = TRUE)
    
    rmse <- sapply(ebnm_res, function(res) {
      return(sqrt(mean((coef(res) - sim_data$theta)^2)))
    })
    
    confint_cov <- sapply(ebnm_res, function(res) {
      zz <- capture.output({ # Capture horseshoe output
        ci <- confint(res, level = 0.9)
      })
      return(1 - mean(sim_data$theta < ci[, 1] | sim_data$theta > ci[, 2]))
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

saveRDS(all_res, "../output/simstudy.rds")

#### Create results table ----

all_res <- all_res %>%
  mutate(Function = str_remove_all(Function, "ebnm_"))

ebnm_fns <- unique(all_res$Function)
sim_fns <- unique(all_res$SimFn)

res_table <- all_res %>%
  select(-SimNumber) %>%
  group_by(Function, SimFn) %>%
  summarize_all(mean) %>%
  ungroup() %>%
  mutate(
    Function = factor(Function, levels = ebnm_fns),
    SimFn = factor(SimFn, levels = sim_fns)
  ) %>%
  arrange(SimFn, Function) %>%
  rename(
    LogLik = LogLikelihood,
    CICov = ConfIntCov
  ) %>%
  pivot_wider(
    names_from = SimFn, 
    values_from = LogLik:CICov,
    names_sep = "XXX"
  )

lbls <- c("Prior Family", rep(c("Normal", "Point-t", "Tophat"), times = 3))
names(lbls) <- names(res_table)
lbls <- as.list(lbls)

rmse_pal <- function(x) {
  min_x <- min(x)
  f <- scales::col_numeric(
    palette = "Reds",
    domain = c(min_x, 1),
    reverse = FALSE
  )
  ifelse(x > 1, f(1), f(x))
}

ci_pal <- function(x) {
  pal_width <- max(abs(range(x) - 0.9))
  f_undercover <- scales::col_numeric(
    palette = "Greens",
    domain = c(0.9 - pal_width, 0.9),
    reverse = TRUE
  )
  f_overcover <- scales::col_numeric(
    palette = "Blues",
    domain = c(0.9, 0.9 + pal_width),
    reverse = FALSE
  )
  ifelse(x < 0.9, f_undercover(x), f_overcover(x))
}

tbl <- res_table %>%
  gt() %>%
  tab_spanner(
    label = "Log likelihood",
    columns = starts_with("Log")
  ) %>%
  tab_spanner(
    label = "RMSE",
    columns = starts_with("RMSE")
  ) %>%
  tab_spanner(
    label = "CI coverage",
    columns = starts_with("CICov")
  ) %>%
  fmt_number(
    columns = starts_with("Log"),
    decimals = 1
  ) %>%
  fmt_number(
    columns = starts_with(c("RMSE", "CI")),
    n_sigfig = 3 
  ) %>%
  data_color(
    columns = starts_with("Log"),
    colors = scales::col_numeric(
      palette = "Reds",
      domain = NULL,
      reverse = TRUE
    )) %>%
  data_color(
    columns = starts_with("RMSE"),
    colors = rmse_pal) %>%
  data_color(
    columns = starts_with("CI"),
    colors = ci_pal
  ) %>%
  cols_label(.list = lbls) %>%
  cols_align("left", columns = Function) %>%
  tab_style(
    style = "padding-left:12px;padding-right:12px;",
    locations = cells_body()
  )

gtsave(tbl, "../figs/simres.htm")
# convert using Adobe: width 9", height 5", margins 0.1", landscape, scale to fit


