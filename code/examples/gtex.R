library(tidyverse)
library(flashier)

url_prefix <- "https://github.com/stephenslab/gtexresults/blob/master/data/"
gtex_url <- paste0(url_prefix, "MatrixEQTLSumStats.Portable.Z.rds?raw=TRUE")
colors_url <- paste0(url_prefix, "GTExColors.txt?raw=TRUE")

gtex <- readRDS(gzcon(url(gtex_url)))
strong <- t(gtex$strong.z) # Dataset used by Urbut et al. and Wang & Stephens.

gtex.colors <- read_tsv(colors_url, col_names = c("Tissue", "Hex", "RGB")) %>%
  mutate(Tissue = str_remove_all(Tissue, "[\\(\\)\\-]")) %>%
  mutate(Tissue = str_replace_all(Tissue, " +", "_")) %>%
  pull(Hex, name = Tissue)
gtex.colors <- gtex.colors[rownames(strong)]

# Point-normal factorization.
pn_res <- flash.init(strong, S = 1, var.type = 2) %>%
  flash.add.greedy(Kmax = 50, ebnm.fn = ebnm::ebnm_point_normal) %>%
  flash.backfit() %>%
  flash.nullcheck()

# Semi-nonnegative factorization.
snn_res <- flash.init(strong, S = 1, var.type = 2) %>%
  flash.add.greedy(
    Kmax = 50, 
    ebnm.fn = c(ebnm::ebnm_point_exponential, ebnm::ebnm_point_normal),
    init.fn = function(f) init.fn.default(f, dim.signs = c(1, 0))) %>%
  flash.backfit() %>%
  flash.nullcheck()

# Extract loadings (tissues).
pn_LL <- pn_res$L.pm[, order(-pn_res$pve)]
snn_LL <- snn_res$L.pm[, order(-snn_res$pve)]

# Flip loadings so that the largest loadings are positive.
pn_sign <- 2L * (apply(pn_LL, 2, max) > -apply(pn_LL, 2, min)) - 1
pn_LL <- t(t(pn_LL) * pn_sign)

# l1-normalize loadings.
pn_LL <- scale(pn_LL, center = FALSE, scale = apply(pn_LL, 2, max))
snn_LL <- scale(snn_LL, center = FALSE, scale = apply(snn_LL, 2, max))

# Rearrange factors to make side-by-side comparison easier.
LL_cor <- abs(cor(pn_LL, snn_LL))
pn_LL <- pn_LL[, apply(LL_cor, 2, which.max)]

colnames(pn_LL) <- paste0("Factor", 1:ncol(pn_LL))
colnames(snn_LL) <- paste0("Factor", 1:ncol(snn_LL))

pn_tib <- as_tibble(pn_LL) %>% add_column(PriorFamily = "point-normal")
snn_tib <- as_tibble(snn_LL) %>% add_column(PriorFamily = "semi-nonnegative")
tib <- pn_tib %>% bind_rows(snn_tib)

tib <- tib %>%
  mutate(Tissue = rep(rownames(strong), 2)) %>%
  pivot_longer(
    -c(PriorFamily, Tissue), names_to = "Factor", values_to = "Loading"
  ) %>%
  mutate(Factor = as.numeric(str_remove(Factor, "Factor")))

plt <- ggplot(tib, aes(x = Tissue, y = Loading, fill = Tissue)) +
  geom_col() +
  facet_grid(rows = vars(Factor), cols = vars(PriorFamily)) +
  scale_fill_manual(values = gtex.colors) +
  scale_y_continuous(breaks = c(0, 0.5, 1)) +
  theme_minimal() + 
  labs(x = "\nTissue", y = "Loading (l1-normalized)\n") +
  guides(fill = guide_legend(ncol = 4)) +
  theme(legend.position = "bottom",
        legend.key.size = unit(0.5, "lines"),
        legend.text = element_text(size = 6.5),
        legend.title = element_blank(),
        axis.text.x = element_blank())

ggsave("../../figs/gtex.png", height = 10, width = 7.5)
