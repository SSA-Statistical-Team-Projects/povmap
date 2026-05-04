
# file tests/testthat/test-megb.R:

library(testthat)

# test data
data(population_data)
data(sample_data)

# tests for megb()
test_that("megb() returns object of class MEGB", {
  result <- megb(
    Y = sample_data$y,
    X = sample_data[, c("x1", "x2", "x3", "x4", "x5")],
    dom_name = "domain",
    smp_data = sample_data,
    pop_data = population_data,
    gradient_params = list(eta = 0.01, nrounds = 50, max_depth = 3),
    seed = 63,
    mse = FALSE
  )
  expect_s3_class(result, "MEGB")
  expect_named(result, c("call", "Indicators", "gbmodel", "inp_smp_data", "unit_preds_all",
                         "unit_pred_smp", "MSE_Estimates", "time_gradient", "gradient_params"))
})

test_that("megb() calculates predictions correctly", {
  result <- megb(
    Y = sample_data$y,
    X = sample_data[, c("x1", "x2", "x3", "x4", "x5")],
    dom_name = "domain",
    smp_data = sample_data,
    pop_data = population_data,
    gradient_params = list(eta = 0.01, nrounds = 50, max_depth = 3),
    seed = 63,
    mse = FALSE
  )
  expect_false(any(is.na(result$unit_pred_smp)))
  expect_false(any(is.na(result$unit_preds_all$unit_preds)))
})

test_that("megb() results in error if sample- and population data are unequal", {
  wrong_population_data <- population_data
  wrong_population_data$domain <- wrong_population_data$domain + 100
  expect_error(
    megb(
      Y = sample_data$y,
      X = sample_data[, c("x1", "x2", "x3", "x4", "x5")],
      dom_name = "domain",
      smp_data = sample_data,
      pop_data = wrong_population_data,
      gradient_params = list(eta = 0.01, nrounds = 50, max_depth = 3),
      seed = 63
    ),
    "There are domains in smp_data that are not present in pop_data."
  )
})
