# test data
# has to be checked

test_that("sample_data and population_data are available and valid", {
  # check if data exists in namespace
  expect_true(exists("sample_data", envir = asNamespace("MEGB")))
  expect_true(exists("population_data", envir = asNamespace("MEGB")))
  
  # get through namespace 
  sample_data <- get("sample_data", envir = asNamespace("MEGB"))
  population_data <- get("population_data", envir = asNamespace("MEGB"))
  
  # structure-tests
  expect_s3_class(sample_data, "data.frame")
  expect_s3_class(population_data, "data.frame")
  
  # check variable names
  expect_true(all(c("x1", "y", "domain") %in% names(sample_data)))
  expect_true(all(c("x1", "y", "domain") %in% names(population_data)))
  
  # no empty data
  expect_gt(nrow(sample_data), 0)
  expect_gt(nrow(population_data), 0)
  
  # check for NAs
  expect_false(anyNA(sample_data))
  expect_false(anyNA(population_data))
})

# is megb reproducible with fixed seed
test_that("megb is reproducible with fixed seed", {
  params <- list(eta = 0.1, nrounds = 100)
  set.seed(1)
  res1 <- megb(
    Y = sample_data$y,
    X = sample_data[, 2:6],
    dom_name = "domain",
    smp_data = sample_data,
    pop_data = population_data,
    gradient_params = params,
    na.rm = TRUE,
    seed = 1,
    mse = TRUE,
    B = 2,
    bootstrap_cores = 0,
    gbm_engine = "xgboost"
  )
  set.seed(1)
  res2 <- megb(
    Y = sample_data$y,
    X = sample_data[, 2:6],
    dom_name = "domain",
    smp_data = sample_data,
    pop_data = population_data,
    gradient_params = params,
    na.rm = TRUE,
    seed = 1,
    mse = TRUE,
    B = 2,
    bootstrap_cores = 0,
    gbm_engine = "xgboost"
  )
  # compare results for one of the output elements:
  if ("unit_preds_all" %in% names(res1)) {
    expect_equal(res1$unit_preds_all, res2$unit_preds_all, tolerance = 1e-12)
  } else {
    expect_equal(res1, res2, tolerance = 1e-12)
  }
})