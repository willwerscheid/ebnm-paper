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
ggsave("./figs/wOBA_comp_cdf.pdf", height = 3, width = 10, units = "in")

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
  subset = most_PA,
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

flash_n <- flash(gtex, ebnm_fn = ebnm_normal, backfit = TRUE)
plot(flash_n, include_scree = FALSE, pm_colors = gtex_colors) +
  ggtitle("") +
  guides(fill = guide_legend(title = "", nrow = 16)) +
  theme(legend.text = element_text(size = 9),
        legend.key.size = unit(6, "points"),
        legend.position = "bottom") 
ggsave("../figs/gtex_n.pdf", height = 7, width = 8, units = "in")

flash_pn <- flash(gtex, ebnm_fn = ebnm_point_normal, backfit = TRUE)
plot(flash_pn, include_scree = FALSE, pm_colors = gtex_colors) +
  ggtitle("") +
  guides(fill = guide_legend(title = "", nrow = 16)) +
  theme(legend.text = element_text(size = 9),
        legend.key.size = unit(6, "points"),
        legend.position = "bottom") 
ggsave("../figs/gtex_pn.pdf", height = 8, width = 8, units = "in")

flash_snn <- flash(gtex, 
                   ebnm_fn = c(ebnm_point_normal, ebnm_point_exponential), 
                   backfit = TRUE)
plot(flash_snn, include_scree = FALSE, pm_colors = gtex_colors) +
  ggtitle("") +
  guides(fill = guide_legend(title = "", nrow = 16)) +
  theme(legend.text = element_text(size = 9),
        legend.key.size = unit(6, "points"),
        legend.position = "bottom") 
ggsave("../figs/gtex_pe.pdf", height = 8, width = 8, units = "in")
