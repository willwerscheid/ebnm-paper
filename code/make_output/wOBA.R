cat("\nCalculating wOBA using FanGraphs data.\n\n")

# Season data is exported from FanGraphs and must include AB, BB, IBB, HBP, 
#   1B, 2B, 3B, HR, and SF. Pitcher data is used to remove pitchers from
#   consideration and need only include IP. wOBA weights are at
#   https://www.fangraphs.com/guts.aspx

calc_wOBA <- function(season, hitter_stats, pitcher_stats = NULL) {
  wts <- read_csv("../../data/wOBA_wts.csv") %>%
    filter(Season == season)
  
  if (!is.null(pitcher_stats)) {
    hitter_stats <- hitter_stats %>%
      left_join(pitcher_stats %>% select(Name, IP), by = "Name") %>%
      filter(is.na(IP) | AB > IP) %>%
      select(-IP)
  }
  
  stats <- hitter_stats %>%
    mutate(uBB = BB - IBB) %>%
    rename(n1B = `1B`, n2B = `2B`, n3B = `3B`) %>%
    mutate(n = uBB + HBP + AB + SF,
           Out = n - uBB - HBP - n1B - n2B - n3B - HR) %>%
    filter(n > 0)
  
  wOBA_wts <- wts %>% 
    mutate(wOut = 0) %>%
    select(wBB, wHBP, w1B, w2B, w3B, wHR, wOut) %>%
    as.matrix()
  
  wOBA_mat <- stats %>% 
    select(uBB, HBP, n1B, n2B, n3B, HR, Out, n) %>%
    mutate_all(~ . / n) %>%
    select(-n) %>%
    as.matrix()
  
  lg_totals <- stats %>% 
    select(uBB, HBP, n1B, n2B, n3B, HR, Out) %>%
    summarize_all(sum) %>%
    as.matrix()
  lg_probs <- lg_totals / sum(lg_totals)
  lg_wOBA <- drop(lg_probs %*% t(wOBA_wts))
  lg_wOBA_sd <- drop(sqrt(lg_probs %*% t(wOBA_wts^2) - lg_wOBA^2))
  
  # wOBA calculated here differs very slightly from published wOBA (but never 
  #   by more than a point), likely because a different precision has been used 
  #   for the linear weights.
  stats <- stats %>%
    mutate(wOBA = drop(wOBA_mat %*% t(wOBA_wts)),
           wOBA_sd = drop(sqrt(wOBA_mat %*% t(wOBA_wts^2) - wOBA^2)) / sqrt(n))
  
  # Use sd based on league probabilities as a floor for the wOBA sds.
  stats <- stats %>%
    mutate(wOBA_sd = pmax(wOBA_sd, lg_wOBA_sd / sqrt(n)))
  
  stats <- stats %>%
    select(Name, playerid, n, wOBA, wOBA_sd)
  
  return(stats)
}

wOBA2020 <- calc_wOBA(2020, 
                      read_csv("../../data/hitters_2020.csv"))
saveRDS(wOBA2020, "../../output/wOBA2020.rds")

for (yr in 2019:2018) {
  wOBA <- calc_wOBA(yr, 
                    read_csv(paste0("../../data/hitters_", yr, ".csv")), 
                    read_csv(paste0("../../data/pitchers_", yr, ".csv")))
  saveRDS(wOBA, paste0("../../output/wOBA", yr, ".rds"))
}
