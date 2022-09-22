## Quality of ashr grid approximations.

# Scale mixtures of normals.

tib <- ebnm:::smngrid
ggplot(tib, aes(x = m, y = KL)) +
  geom_line() +
  scale_y_log10() +
  labs(x = "\nGrid multiplier", y = "KL divergence (upper bound)\n") +
  theme_minimal()
ggsave("../../figs/smnKLdiv.png", width = 6, height = 3.375)

# Symmetric unimodals.

KLs <- 10^c(-2, -2.5, -3, -4)

symmunigrid <- ebnm:::symmunigrid %>%
  filter(KL %in% KLs)

tib <- symmunigrid %>%
  left_join(symmunigrid %>% mutate(idx = idx - 1), by = c("KL", "idx"), suffix = c("", ".next")) %>%
  mutate(grid.ratio.m1 = loc.next / loc - 1) %>%
  mutate(theory.lim = exp(1) * KL) %>%
  mutate(KL = signif(KL, digits = 1)) %>%
  mutate(KL = factor(paste("Target upper bound:", KL),
                     levels = paste("Target upper bound:", -sort(-unique(KL))))) %>%
  filter(idx > 1, !is.na(loc.next))

ggplot(tib, aes(x = idx, y = grid.ratio.m1)) +
  geom_point(size = 0.25) +
  geom_line(aes(x = idx, y = theory.lim), linetype = "dashed") +
  facet_wrap(~KL) +
  scale_y_log10() +
  labs(x = "\nGrid index (k)", y = TeX("Grid ratio $\\frac{a_k}{a_{k-1}}$ - 1")) +
  theme_minimal()
ggsave("../../figs/symmKLdiv.png", width = 6, height = 3.375)

# KL divergence for default ashr grid (symmetric unimodal priors).

x.max <- 10^seq(log10(3), log10(50), length.out = 40)
ub <- sapply(x.max, function(x) {
  ebnm:::ub.symmKLdiv(x / sqrt(2), x)$KL.div
})
tib <- tibble(x.max = x.max, ub = ub)

ggplot(tib, aes(x = x.max, y = ub)) +
  geom_line() +
  scale_y_log10() +
  labs(x = "\nMaximum relevant grid point", y = "KL divergence (upper bound)\n") +
  theme_minimal()
ggsave("../../figs/ashr_symmKLdiv.png", width = 6, height = 3.375)

# NPMLE.

d_over_s <- 10^seq(-1, 1, length.out = 50)
tib <- tibble()
for (d in d_over_s) {
  optres <- ebnm:::min.npmleKLdiv(d)
  tib <- tib %>%
    bind_rows(tibble(d_over_s = d, tight = optres$objective, s2 = optres$minimum))
}

ez.bound <- function(x) pmin(x^4 / 64, log(1 + x^2 / 4) / 2)
tib <- tib %>% mutate(simple = ez.bound(d_over_s))
tib <- tib %>%
  pivot_longer(c(tight, simple), names_to = "ub.type", values_to = "ub")

ggplot(tib, aes(x = d_over_s, y = ub)) +
  geom_line(aes(linetype = ub.type)) +
  scale_y_log10() +
  labs(x = "\nDistance between grid points (scaled by SE)",
       y = "KL divergence (upper bound)\n",
       linetype = "Bound type") +
  theme_minimal() +
  scale_linetype_manual(values = c("dashed", "solid"))
ggsave("../../figs/npmleKLdiv.png", width = 6, height = 3.375)
