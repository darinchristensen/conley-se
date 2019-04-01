context("test-conley")
library("ConleySE")
library("lfe")

# Loading test data:
data("conley_spatial")

test_that("OLS model using lm() works", {
  cs2012 <- conley_spatial[conley_spatial$year == 2012,]
  lm1 <- lm(EmpClean00 ~ HDD + CDD, data = cs2012)
  
  felm1 <- felm(EmpClean00 ~ HDD + CDD | 0 | 0 | latitude + longitude,
            data = cs2012, keepCX = TRUE)
  
  
  expected_res <- c(56599.358, 12.971, 35.198)
  names(expected_res) <- c("(Intercept)", "HDD", "CDD")
  lm_vcov_spatial <- ConleySE(model = lm1,
                 x = cs2012$longitude,
                 y = cs2012$latitude,
                 dist_fn = "SH", dist_cutoff = 500, 
                 verbose = FALSE) 
  
  felm_vcov_spatial <- ConleySE(model = felm1,
                                x_var = "longitude",
                                y_var = "latitude",
                              dist_fn = "SH", dist_cutoff = 500, 
                              verbose = FALSE) 
  
  
  lm_se_spatial <- round(sqrt(diag(lm_vcov_spatial)), 3)
  felm_se_spatial <- round(sqrt(diag(felm_vcov_spatial$Spatial)), 3)
  expect_equal(se_spatial, expected_res)
  expect_equal(felm_se_spatial, expected_res)
})

test_that("Fixed effects model using plm::plm() matches expected results", {
  if(require(plm)) {
    
    pdf <- pdata.frame(conley_spatial, index = c("FIPS", "year"), drop.index = F, row.names = T)
    pm1 <- plm(EmpClean00 ~ HDD + CDD, data = pdf[pdf$year == 2012,], model = "pooling")
    
    G <- length(unique(pdf$FIPS))
    N <- nrow(pdf)
    dfa <- (G/(G - 1)) * (N - 1)/pm1$df.residual
    
    # display with cluster VCE and df-adjustment
    firm_c_vcov <- dfa * vcovHC(pm1, type = "HC0", cluster = "group", adjust = T)
    coeftest(pm1, vcov = firm_c_vcov)
    
  } else {
    message("Error: package 'plm' not installed.")
    FALSE
  }
})

test_that("Fixed effects model using lfe::felm() matches expected results", {
  if(require(lfe)) {
    # Note that the original test results are for mislabelled latitude and
    # longitude variables, so we have to rename them here for the test to 
    # work and match results from Darin Christensen's repo. The data in this 
    # package are correctly labelled.
    names(conley_spatial)[6:7] <- c("latitude", "longitude")
    
    expected_res <- matrix(c(0.650, 1.493, 0.886, 4.065, 0.721, 3.631), 
                           nrow = 2,
                           byrow = FALSE)
    colnames(expected_res) <- c("OLS", "Spatial", "Spatial_HAC")
    rownames(expected_res) <- c("HDD", "CDD")
    
    m <- felm(EmpClean00 ~ HDD + CDD | year + FIPS | 0 | latitude + longitude,
              data = conley_spatial[!is.na(conley_spatial$EmpClean00),], keepCX = TRUE)
    
    SE <- ConleySE(model = m,
                   id_var = "FIPS", 
                   time_var = "year",
                   x_var = "longitude", y_var = "latitude",
                   dist_fn = "SH", dist_cutoff = 500, 
                   lag_cutoff = 5,
                   cores = 1, 
                   verbose = FALSE) 
    res <- sapply(SE, function(x) round(diag(sqrt(x)), 3))
    expect_equal(res, expected_res)
  } else {
    message("Error: package 'lfe' not installed.")
  }
})
