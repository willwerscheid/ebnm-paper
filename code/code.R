###### HANDLE ARGS ----------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

valid.args <- args %in% c(
  "out-to-file", 
  "test",
  "full"
)
if (!all(valid.args)) {
  stop("Command line argument ",
       min(which(!valid.args)),
       " not recognized.")
}

test <- "test" %in% args
full <- "full" %in% args

if ("out-to-file" %in% args) {
  fname <- "../output/code_output"
  if (test) {
    fname <- paste0(fname, "_test")
  } 
  out_file <- file(paste0(fname, ".txt"), open = "wt")
  sink(out_file)
  sink(out_file, type = "message")
}

system("if [ ! -d ../figs ]; then mkdir ../figs; fi")


###### REQUIRED PACKAGES ----------------------------------------------

library(tidyverse)
library(ebnm)
library(flashier)
library(gt)

# Also needed: microbenchmark, scales, Rtsne, ggrepel, cowplot


###### TIMING STUDY ---------------------------------------------------

cat("\nTIMING STUDY...\n\n")

#### Run simulations ----

set.seed(666)

if (exists("test") && test) {
  ns <- 10^seq(2, 3, by = 0.5)
  mb_times <- rep(2, length(ns))
  fname_suffix <- "_test"
} else if (exists("full") && full) {
  ns <- 10^seq(2, 6, by = 0.5)
  mb_times <- rep(20, length(ns))
  fname_suffix <- "_full"
} else {
  ns <- 10^seq(2, 4, by = 0.5)
  mb_times <- rep(20, length(ns))
  fname_suffix <- ""
}

cat(paste0("\nTiming different prior families.\n",
           "  n ranges from ", min(ns), " to ", max(ns), ".\n\n"))

# Mixture of point mass at zero and t_5
point_t <- function(n) {
  samp <- 1.5 * rt(n, df = 5)
  samp <- samp * sample(c(0, 1), n, replace = TRUE, prob = c(0.5, 0.5))
  return(samp + rnorm(n))
}

ebnm_fns <- c("ebnm_normal",
              "ebnm_point_normal",
              "ebnm_point_laplace",
              "ebnm_normal_scale_mixture",
              "ebnm_unimodal_symmetric",
              "ebnm_unimodal",
              "ebnm_deconvolver",
              "ebnm_npmle",
              "ebnm_horseshoe")

t_begin <- Sys.time()
res <- tibble()
for (i in 1:length(ns)) {
  n <- ns[i]
  cat("  n:", n, "\n")
  
  x <- point_t(n)
  
  test_fns <- ebnm_fns
  
  mb_tests <- lapply(test_fns, function(fn) {
    bquote(do.call(.(fn), list(x = x, s = 1)))
  })
  
  mb_res <- microbenchmark::microbenchmark(
    list = mb_tests,
    times = mb_times[i]
  )
  
  res <- res %>%
    bind_rows(
      tibble(mb_res) %>%
        mutate(expr = str_extract(expr, "ebnm.+\\\"")) %>%
        mutate(expr = str_remove(expr, "\\\"")) %>%
        mutate(n = n)
    )
}

t_elapsed <- Sys.time() - t_begin
cat("  Done. Time elapsed:",
    round(as.numeric(t_elapsed, units = "mins"), 1), "minutes.\n")

res <- res %>%
  mutate(expr = fct_relevel(expr, rev(ebnm_fns)))

saveRDS(res, paste0("../output/timecomps", fname_suffix, ".rds"))

#### Create figure ----

summary_res <- res %>%
  mutate(time = time / 1e9) %>%
  group_by(n, expr) %>%
  summarize(
    mean = mean(time),
    lowerq = quantile(time, probs = 0.1),
    upperq = quantile(time, probs = 0.9)
  ) 

lvls <- summary_res %>%
  filter(n == ns[8]) %>%
  arrange(-mean) %>%
  pull(expr) %>%
  as.character()

summary_res <- summary_res %>%
  mutate(expr = fct_relevel(expr, lvls)) %>%
  mutate(expr = fct_relabel(expr, ~ str_remove(., "ebnm_")))

ggplot(summary_res, aes(x = n, y = mean, color = expr)) +
  geom_point() +
  geom_line(linewidth = 0.2) +
  geom_errorbar(aes(ymin = lowerq, ymax = upperq), width = 0.05) +
  scale_x_log10() +
  scale_y_log10() +
  scale_color_brewer(palette = "Set1") +
  labs(x = "\nNumber of observations", 
       y = "Runtime (s)\n", 
       col = "Prior family") +
  theme_minimal()

ggsave(paste0("../figs/timecomps", fname_suffix, ".pdf"), 
       width = 6, height = 3.375)


###### SIMULATION STUDY -----------------------------------------------

cat("\nSIMULATION STUDY...\n\n")

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
  fname_suffix <- "_test"
} else {
  nsim <- 10
  n <- 1000
  fname_suffix <- ""
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

saveRDS(all_res, paste0("../output/simstudy", fname_suffix, ".rds"))

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

gtsave(tbl, paste0("../figs/simres", fname_suffix, ".htm"))
# convert using Adobe: width 9", height 5", margins 0.1", landscape, scale to fit


###### EXAMPLES -------------------------------------------------------

cat("\nEXAMPLES...\n\n")

#### wOBA (part I) ----

library(ebnm)
data(wOBA)
nrow(wOBA)
head(wOBA)

library(ggplot2)
ggplot(wOBA, aes(x = x)) +
  geom_histogram(bins = 64, color = "black") +
  theme_classic()
ggsave("../figs/wOBA_dist.pdf", height = 3, width = 5, units = "in")

x <- wOBA$x
s <- wOBA$s
names(x) <- wOBA$Name
names(s) <- wOBA$Name
fit_normal <- ebnm(x, s, prior_family = "normal", mode = "estimate")

fit_normal <- ebnm_normal(x, s, mode = "estimate")

summary(fit_normal)

plot(fit_normal)
ggsave("../figs/wOBA_normal.pdf", height = 3, width = 5, units = "in")

plot(fit_normal) +
  geom_point(aes(color = sqrt(wOBA$PA))) +
  labs(x = "wOBA", y = "EB estimate of true wOBA skill", 
       color = expression(sqrt(PA))) +
  scale_color_gradient(low = "blue", high = "red")
ggsave("../figs/wOBA_normal_custom.pdf", height = 3, width = 5, units = "in")

print(head(fitted(fit_normal)), digits = 3)

fit_unimodal <- ebnm(x, s, prior_family = "unimodal", mode = "estimate")

top50 <- order(wOBA$PA, decreasing = TRUE)
top50 <- top50[1:50]
plot(fit_normal, fit_unimodal, subset = top50)
ggsave("../figs/wOBA_comp.pdf", height = 3, width = 5, units = "in")

dat <- cbind(wOBA[, c("PA","x")],
             fitted(fit_normal),
             fitted(fit_unimodal))
names(dat) <- c("PA", "x", "mean1", "sd1", "mean2", "sd2")
print(head(dat), digits = 3)

library(cowplot)
p1 <- plot(fit_normal, fit_unimodal, incl_cdf = TRUE, incl_pm = FALSE) +
  xlim(c(.250, .350)) +
  guides(color = "none")
p2 <- plot(fit_normal, fit_unimodal, incl_cdf = TRUE, incl_pm = FALSE) +
  lims(x = c(.350, .450), y = c(0.95, 1))
plot_grid(p1, p2, nrow = 1, ncol = 2, rel_widths = c(.4, .6))
ggsave("../figs/wOBA_comp_cdf.pdf", height = 3, width = 10, units = "in")

fit_unimodal <- ebnm_add_sampler(fit_unimodal)
set.seed(1)
print(head(confint(fit_unimodal, level = 0.8)), digits = 3)


#### wOBA (part II) ----

fit_npmle <- ebnm(x, s, prior_family = "npmle")

fit_npmle <- ebnm(x, s, prior_family = "npmle", 
                  control = list(verbose = TRUE))

plot(fit_normal, fit_unimodal, fit_npmle, incl_cdf = TRUE, subset = top50)

# Slightly different from the text (need to show one plot at a time):
plot(
  fit_normal, fit_unimodal, fit_npmle,
  incl_cdf = TRUE,
  incl_pm = FALSE
) + xlim(0.25, 0.45) 
ggsave("../figs/wOBA_npmle_cdf.pdf", height = 4, width = 6, units = "in")
plot(
  fit_normal, fit_unimodal, fit_npmle,
  subset = top50,
  incl_cdf = FALSE,
  incl_pm = TRUE
)
ggsave("../figs/wOBA_npmle_pm.pdf", height = 4, width = 6, units = "in")

logLik(fit_unimodal)
logLik(fit_npmle)

scale_npmle <- ebnm_scale_npmle(x, s, KLdiv_target = 0.001/length(x), 
                                max_K = 1000)
fit_npmle_finer <- ebnm_npmle(x, s, scale = scale_npmle)
logLik(fit_npmle)
logLik(fit_npmle_finer)

fit_npmle <- ebnm_add_sampler(fit_npmle)
print(head(quantile(fit_npmle, probs = c(0.1, 0.9))), digits = 3)

confint(fit_npmle, level = 0.8, parm = "Aaron Judge")

fit_deconv <- ebnm_deconvolver(x / s, output = ebnm_output_all()) 
plot(fit_deconv, incl_cdf = TRUE, incl_pm = FALSE)
ggsave("../figs/wOBA_deconv.pdf", height = 4, width = 6, units = "in")

print(head(quantile(fit_deconv, probs = c(0.1, 0.9)) * s), digits = 3)


#### GTEx ----

library(flashier)
data(gtex)
nrow(gtex)
ncol(gtex)
gtex[1:2, 1:2]

library(Rtsne)
library(ggrepel)
set.seed(1)
out <- Rtsne(t(gtex), dims = 2, perplexity = 10)
pdat <- data.frame(d1 = out$Y[, 1],
                   d2 = out$Y[, 2],
                   tissue = colnames(gtex))
ggplot(pdat,aes(x = d1, y = d2, label = tissue)) +
  geom_point(size = 6, color = gtex_colors) +
  geom_text_repel(size = 2.5, max.overlaps = Inf) +
  labs(x = "t-SNE [1]", y = "t-SNE [2]") +
  theme_classic() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank())
ggsave("../figs/gtex_tsne.pdf", height = 4, width = 6, units = "in")

t_n <- system.time({
  flash_n <- flash(gtex, ebnm_fn = ebnm_normal, backfit = TRUE)
})
t_n[3]

plot(flash_n, include_scree = FALSE, pm_colors = gtex_colors) +
  ggtitle("") +
  guides(fill = guide_legend(title = "", nrow = 16)) +
  theme(legend.text = element_text(size = 9),
        legend.key.size = unit(6, "points"),
        legend.position = "bottom") 
ggsave("../figs/gtex_n.pdf", height = 7, width = 8, units = "in")

t_pn <- system.time({
  flash_pn <- flash(gtex, ebnm_fn = ebnm_point_normal, backfit = TRUE)
})
t_pn[3]

plot(flash_pn, include_scree = FALSE, pm_colors = gtex_colors) +
  ggtitle("") +
  guides(fill = guide_legend(title = "", nrow = 16)) +
  theme(legend.text = element_text(size = 9),
        legend.key.size = unit(6, "points"),
        legend.position = "bottom") 
ggsave("../figs/gtex_pn.pdf", height = 8, width = 8, units = "in")

t_snn <- system.time({
  flash_snn <- flash(gtex, 
                     ebnm_fn = c(ebnm_point_normal, ebnm_point_exponential), 
                     backfit = TRUE)
})
t_snn[3]

plot(flash_snn, include_scree = FALSE, pm_colors = gtex_colors) +
  ggtitle("") +
  guides(fill = guide_legend(title = "", nrow = 16)) +
  theme(legend.text = element_text(size = 9),
        legend.key.size = unit(6, "points"),
        legend.position = "bottom") 
ggsave("../figs/gtex_pe.pdf", height = 8, width = 8, units = "in")


###### SESSION INFO ---------------------------------------------------

cat("\n\n")
sessionInfo()
