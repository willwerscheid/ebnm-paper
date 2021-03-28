## wOBA example.

cat("\nwOBA example...\n\n")

wOBA2020 <- readRDS("../../output/wOBA2020.rds")
wOBA2019 <- readRDS("../../output/wOBA2019.rds")
wOBA2018 <- readRDS("../../output/wOBA2018.rds")

wOBA_threeyr <- wOBA2020 %>% 
  full_join(wOBA2019, by = "playerid", suffix = c("", "19")) %>% 
  full_join(wOBA2018, by = "playerid", suffix = c("20", "18")) %>%
  mutate(Name = case_when(!is.na(Name20) ~ Name20,
                          !is.na(Name19) ~ Name19,
                          TRUE ~ Name18)) %>%
  mutate_at(vars(wOBA20, wOBA19, wOBA18), replace_na, 0) %>%
  mutate_at(vars(wOBA_sd20, wOBA_sd19, wOBA_sd18), replace_na, Inf) %>%
  mutate(wOBA_sd19 = wOBA_sd19 * 5 / 4,
         wOBA_sd18 = wOBA_sd18 * 5 / 3) %>%
  mutate(wOBA_sd = sqrt(1 / (1 / wOBA_sd20^2 + 1 / wOBA_sd19^2 + 1 / wOBA_sd18^2))) %>%
  mutate(wOBA = (wOBA20 / wOBA_sd20^2 + wOBA19 / wOBA_sd19^2 + wOBA18 / wOBA_sd18^2) * wOBA_sd^2) %>%
  select(Name, playerid, wOBA, wOBA_sd)


# Estimate posterior means using single seasons' worth of data.

ebnm2020 <- ebnm_unimodal(wOBA2020$wOBA, wOBA2020$wOBA_sd, 
                          mode = "estimate", 
                          output = output_all())
wOBA2020 <- wOBA2020 %>% 
  add_column(year = "2020", postmean = ebnm2020$posterior$mean)

ebnm2019 <- ebnm_unimodal(wOBA2019$wOBA, wOBA2019$wOBA_sd, 
                          mode = "estimate", 
                          output = output_all())
wOBA2019 <- wOBA2019 %>% 
  add_column(year = "2019", postmean = ebnm2019$posterior$mean)

ggplot(wOBA2020 %>% bind_rows(wOBA2019), 
       aes(x = wOBA, y = postmean, col = n)) +
  geom_point() +
  geom_abline(slope = 1, linetype = "dashed") +
  labs(x = "wOBA (observed)", y = "wOBA skill (estimated)", col = "PA") +
  facet_wrap(~year, nrow = 2) +
  theme_minimal()

ggsave("../../figs/wOBAbyyear.png", width = 6, height = 3.375)


# Sample from the posterior to estimate hitters' quality.

set.seed(1)
samp2020 <- ebnm2020$posterior_sampler(5000)
samp2019 <- ebnm2019$posterior_sampler(5000)

hitter_quality <- function(samp) {
  c("excellent" = mean(samp > .4),
    "great" = mean(samp > .37),
    "above_avg" = mean(samp > .34))
}

wOBA2020 <- wOBA2020 %>% 
  bind_cols(as_tibble(t(apply(samp2020, 2, hitter_quality))))
wOBA2019 <- wOBA2019 %>% 
  bind_cols(as_tibble(t(apply(samp2019, 2, hitter_quality))))

wOBA2020 %>% 
  bind_rows(wOBA2019) %>% 
  group_by(year) %>% 
  summarize_at(vars(excellent:above_avg), ~sum(. > .8)) %>%
  as.data.frame()


# Use three years' worth of data to estimate current ability.
cat("Generating estimates of current ability:\n\n")

ebnm_threeyr <- ebnm_unimodal(wOBA_threeyr$wOBA, wOBA_threeyr$wOBA_sd, 
                              mode = "estimate", 
                              output = output_all())
wOBA_threeyr <- wOBA_threeyr %>% 
  add_column(postmean = ebnm_threeyr$posterior$mean)

samp_threeyr <- ebnm_threeyr$posterior_sampler(5000)
wOBA_threeyr <- wOBA_threeyr %>% 
  bind_cols(as_tibble(t(apply(samp_threeyr, 2, hitter_quality))))

wOBA_threeyr %>% summarize_at(vars(excellent:above_avg), ~sum(. > .8))

wOBA_threeyr %>% 
  filter(above_avg > .8) %>% 
  select(Name, wOBA, wOBA_sd, postmean:above_avg) %>%
  rename(p_excellent = excellent, p_great = great, p_above_avg = above_avg) %>%
  arrange(-p_above_avg) %>%
  xtable(digits = 3) %>%
  print(include.rownames = FALSE)
