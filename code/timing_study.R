###### TIMING STUDY ---------------------------------------------------

cat("\nTIMING STUDY...\n\n")

#### Run simulations ----

set.seed(666)

if (exists("test") && test) {
  ns <- 10^seq(2, 3, by = 0.5)
  mb_times <- rep(2, length(ns))
} else if (exists("full") && full) {
  ns <- 10^seq(2, 6, by = 0.5)
  mb_times <- rep(20, length(ns))
} else {
  ns <- 10^seq(2, 4, by = 0.5)
  mb_times <- rep(20, length(ns))
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
    cat("    fn:", fn, "\n")
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

saveRDS(res, "./output/timecomps.rds")

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
  labs(x = "\nNumber of observations", y = "Time elapsed (s)\n", col = "Prior Family") +
  theme_minimal()

ggsave("./figs/timecomps.pdf", width = 6, height = 3.375)
