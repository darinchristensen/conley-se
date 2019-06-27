context("test-cpp")
library(ConleySE)

data("conley_spatial")

test_that("DistMatrix() works with 'SH' distance function", {
  x <- c(-120,-115, -112)
  y <- c(45, 44.5, 42.1)
  a <- cbind(x, y)
  d <- DistMatrix(a, cutoff = 500, dist_fn = "SH")
  #print(d)
  res <- matrix(c(1,.2073015,0,.2073015,1,.2861902,0,.2861902,1), nrow = 3)
  expect_equal(round(d,5), round(res,5))
})

test_that("DistMatrix() works with 'Haversine' distance function", {
  cutoff <- 1000
  x <- c(-120, -115, -112)
  y <- c(45, 44.5, 42.1)
  a <- cbind(x, y)
  distance_true <- geosphere::distHaversine(a, rbind(a[-1,], a[1,]))/1000 # in km 
  distance_rel <- (1 - distance_true / cutoff) * (distance_true <= cutoff) # relative distance ("Bartlett")
  d <- DistMatrix(a, cutoff = 1000, dist_fn = "Haversine")
  d_vec <- d[upper.tri(d)]
  expect_equal(round(sort(d_vec), 5), round(sort(distance_rel),5))
})

test_that("Bal_XeeXhC()", {
})

test_that("XeeXhC()", {
})

test_that("XeeXhC_Lg()", {
})

test_that("TimeDist()", {
})



















