# em_maxit / em_tol: the EM cap and stopping tolerance are arguments (defaults 25 and 1e-04, the earlier constants).
test_that("the defaults are the earlier constants", {
  expect_identical(formals(megb)$em_maxit, povmap:::.megb_em_max_iterations)
  expect_identical(formals(megb)$em_tol, povmap:::.megb_em_error_tolerance)
})

n_it <- function(m) m$megb_model$iterations_used

test_that("em_maxit sets the number of EM iterations (point fit and with the bootstrap)", {
  cfg <- megb_bench_cfgs$log_ratio
  for (k in c(2L, 4L)) {
    m <- megb_bench_run(cfg, mse = FALSE, benchmark = NULL, benchmark_level = NULL, em_maxit = k)
    expect_identical(as.integer(n_it(m)), k)
  }
  m <- megb_bench_run(cfg, B = 3, em_maxit = 3L)
  expect_identical(as.integer(n_it(m)), 3L)
})

test_that("em_tol stops the EM early when the relative log-likelihood change is below it", {
  cfg <- megb_bench_cfgs$log_ratio
  m <- megb_bench_run(cfg, mse = FALSE, benchmark = NULL, benchmark_level = NULL, em_maxit = 12L, em_tol = 1e6)
  expect_lt(as.integer(n_it(m)), 12L)
})

test_that("em_maxit and em_tol are validated", {
  cfg <- megb_bench_cfgs$log_ratio
  for (bad in list(0, -1, 2.5, NA, "a", c(5, 6))) expect_error(megb_bench_run(cfg, em_maxit = bad), "em_maxit must be a whole number")
  for (bad in list(0, -1, NA, "a", c(1e-4, 1e-3))) expect_error(megb_bench_run(cfg, em_tol = bad), "em_tol must be a positive number")
})
