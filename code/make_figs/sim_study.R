all_res <- readRDS("../../output/simstudy.rds")

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

lbls <- c("Function", rep(c("LogLik", "RMSE", "CICov"), each = 3))
names(lbls) <- names(res_table)
lbls <- as.list(lbls)

tbl <- res_table %>%
  gt() %>%
  tab_spanner(
    label = "Point-normal",
    columns = ends_with("normal")
  ) %>%
  tab_spanner(
    label = "Point-t",
    columns = ends_with("_t")
  ) %>%
  tab_spanner(
    label = "Asymmetric tophat",
    columns = ends_with("tophat")
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
    colors = scales::col_numeric(
      palette = "Reds",
      domain = NULL,
      reverse = FALSE
    )) %>%
  data_color(
    columns = starts_with("CI"),
    colors = scales::col_numeric(
      palette = "Reds",
      domain = NULL,
      reverse = TRUE
    )) %>%
  cols_label(.list = lbls) %>%
  cols_align("left", columns = Function) %>%
  tab_style(
    style = "padding-left:12px;padding-right:12px;",
    locations = cells_body()
  )

gtsave(tbl, "../../figs/simres.png", vwidth = 900, vheight = 600)
