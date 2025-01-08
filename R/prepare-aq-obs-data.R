library(saqgetr)
library(lubridate)
library(dplyr)
library(readr)
library(glue)

# get annual means, keep one year PM2.5
year <- 2023
adat <- get_saq_simple_summaries(summary = "annual_means")
pm25_eu <- adat %>% filter(variable=="pm2.5", date==ymd(glue("{year}-01-01")))

# get sites infos, merge with data
anag <- get_saq_sites()
left_join(
  pm25_eu,
  anag %>% select(site:site_area,eoi_code)
) -> dat

# keep only background stations with at least 80% of the data
dat %>%
  filter(
    site_type=="background",
    (summary_source==1 & count>=8760*0.8) |   # hourly data
      (summary_source==20 & count>=365*0.8)   # daily data
  ) -> dat

# write in CSV format
write_csv(dat, glue("data/eea_obs_pm25_bkg_{year}.csv"))
