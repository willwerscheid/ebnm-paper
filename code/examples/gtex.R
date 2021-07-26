cat("\nReading GTEx files...\n\n")

gtex <- readRDS(gzcon(url("https://github.com/stephenslab/gtexresults/blob/master/data/MatrixEQTLSumStats.Portable.Z.rds?raw=TRUE")))
strong <- t(gtex$strong.z)

missing.tissues <- c(7, 8, 19, 20, 24, 25, 31, 34, 37)
gtex.colors <- read_tsv("https://github.com/stephenslab/gtexresults/blob/master/data/GTExColors.txt?raw=TRUE", col_names = FALSE, )
gtex.colors <- gtex.colors[-missing.tissues, ]
gtex.colors <- gtex.colors %>% pull(X2)
names(gtex.colors) <- rownames(strong)

cat("\nRunning flash...\n\n")

if (exists("test") && test) {
  greedy.Kmax <- 5
  backfit <- FALSE
} else {
  greedy.Kmax <- 50
  backfit <- TRUE
}

pn_res <- flash(
  strong,
  S = 1,
  greedy.Kmax = greedy.Kmax,
  prior.family = prior.point.normal(),
  var.type = 2,
  backfit = backfit
)

snn_res <- flash(
  strong,
  S = 1,
  greedy.Kmax = greedy.Kmax,
  prior.family = c(prior.nonnegative(), prior.point.normal()),
  var.type = 2,
  backfit = backfit
)

pn_LL <- pn_res$loadings.pm[[1]][, order(-pn_res$pve)]
snn_LL <- snn_res$loadings.pm[[1]][, order(-snn_res$pve)]

# Flip loadings where necessary.
pn_sign <- 2L * (apply(pn_LL, 2, max) > -apply(pn_LL, 2, min)) - 1
pn_LL <- t(t(pn_LL) * pn_sign)

# Rearrange to make the factors line up:
LL_cor <- abs(cor(pn_LL, snn_LL))
snn_LL <- snn_LL[, apply(LL_cor, 1, which.max)]

colnames(pn_LL) <- paste0("Factor", 1:ncol(pn_LL))
colnames(snn_LL) <- paste0("Factor", 1:ncol(snn_LL))
tib <- as_tibble(pn_LL) %>% add_column(PriorFamily = "point-normal") %>%
  bind_rows(as_tibble(snn_LL) %>% add_column(PriorFamily = "semi-nonnegative"))

tib <- tib %>%
  mutate(Tissue = rep(rownames(strong), 2)) %>%
  pivot_longer(-c(PriorFamily, Tissue), names_to = "Factor", values_to = "Loading") %>%
  mutate(Factor = as.numeric(str_remove(Factor, "Factor")))

plot_kset <- function(kset) {
  plt <- ggplot(tib %>% filter(Factor %in% kset), aes(x = Tissue, y = Loading, fill = Tissue)) +
    geom_col() +
    facet_grid(rows = vars(Factor), cols = vars(PriorFamily)) +
    scale_fill_manual(values = gtex.colors) +
    theme_minimal() + 
    labs(x = "", y = "") +
    theme(legend.position = "none",
          axis.text.x = element_blank())
  ggsave(paste0("../../figs/gtex", min(kset), ".png"), height = 6, width = 6)
}

if (exists("test") && test) {
  plot_kset(1:5)
} else {
  plot_kset(1:6)
  plot_kset(7:12)
  plot_kset(13:18)
}

plt <- ggplot(tib %>% filter(PriorFamily == "point-normal",
                             Factor == 1), 
              aes(x = Tissue, y = Loading, fill = Tissue)) +
  geom_col() +
  scale_fill_manual(values = gtex.colors) +
  theme_minimal() + 
  labs(x = "", y = "") +
  guides(fill = guide_legend(ncol = 3)) +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 6.5),
        legend.title = element_blank(),
        axis.text.x = element_blank())

my_legend <- get_legend(plt)
legend_plt <- as_ggplot(my_legend)
ggsave("../../figs/gtex_legend.png", height = 4, width = 7)
