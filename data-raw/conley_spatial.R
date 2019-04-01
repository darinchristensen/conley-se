#!/usr/bin/Rscript
#
# Get and process the source data for testing. Data sourced from the 
# Darin Christensen conley-se repo.
#
# NOTE: For some reason the latitude/longitude variables are misnamed, and so
#       are fixed here. We rename them in the test results so that they match
#       results from Darin Christensen but otherwise leave them corrected in the 
#       data file.

library(dplyr)
fileurl <- "https://github.com/darinchristensen/conley-se/raw/master/data/new_testspatial.dta"

conley_spatial <- foreign::read.dta(fileurl) %>%
  select(FIPS, year, EmpClean00, HDD, CDD, latitude, longitude) %>%
  rename(latitude = longitude, longitude = latitude) %>%
  na.omit()

usethis::use_data(conley_spatial, overwrite = TRUE)
