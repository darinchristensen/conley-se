#!/usr/bin/Rscript
#
# Data on climate and rents in the Western fruitful rim
library(dplyr)

census_vars <- c("x_pop_density", "x_housing_density")
soil_vars <- c("x_salinity", "x_frac_flood", "x_wetland", "x_k_fact", 
               "x_slope_len", "x_sand", "x_clay", "x_moisture", "x_permeability")

frr_climate <- readRDS("data-raw/analysis_linear_0310_huc2_gridmet_cs_ddays.rds") %>%
  filter(in_frr == TRUE) %>%
  filter(sw_irr > 0 | sw_irr_crop > 0) %>%
  filter(wu_peracre < 10) %>%
  mutate(
    ln_y_rent_irr = log(rent_irr),
    gdd = (dday_31 - dday_7)/1000,
    hdd = (dday_40 - dday_34)/1000
  ) %>%
  mutate(gdd_sq = gdd^2,
         hdd_sqrt = sqrt(hdd),
         precip_sq = precip^2) %>%
  select(fips, st_fips, lat, lon, ln_y_rent_irr, 
         gdd, gdd_sq, hdd_sqrt, precip, precip_sq, 
         wu_peracre, acres_crop_harv_irr,
         {{ census_vars }}, {{ soil_vars }}) %>%
  na.omit()

usethis::use_data(frr_climate, overwrite = TRUE)
