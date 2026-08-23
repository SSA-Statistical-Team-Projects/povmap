# Tests for xgb_tune(search = "random") and the early-stopping option.
#
# The governing requirement is that an existing call, with none of the new
# arguments supplied, behaves exactly as it did before random search was added.

make_toy <- function(n_dom = 30, per = 20, seed = 1) {
  set.seed(seed)
  muni <- rep(seq_len(n_dom), each = per)
  x1 <- rnorm(n_dom * per); x2 <- rnorm(n_dom * per); x3 <- rnorm(n_dom * per)
  u  <- rnorm(n_dom, sd = 0.5)[muni]
  data.frame(y = 1 + 0.6 * x1 - 0.4 * x2 + 0.2 * x3 + u + rnorm(n_dom * per, sd = 0.5),
             x1 = x1, x2 = x2, x3 = x3, muni = muni, wt = 1)
}

## ---- helpers: no fitting, no network -------------------------------------

test_that("bare c(lo, hi) inherits the parameter's default scale", {
  expect_equal(.xgb_norm_bound("eta", c(0.01, 0.3))$scale, "log")
  expect_equal(.xgb_norm_bound("min_child_weight", c(1, 50))$scale, "log")
  expect_equal(.xgb_norm_bound("lambda", c(0.1, 10))$scale, "log")
  expect_equal(.xgb_norm_bound("max_depth", c(2, 10))$scale, "int")
  expect_equal(.xgb_norm_bound("alpha", c(0, 5))$scale, "zeroinf")
  expect_equal(.xgb_norm_bound("subsample", c(0.5, 1))$scale, "linear")
})

test_that("the list form overrides the default scale", {
  expect_equal(.xgb_norm_bound("eta", list(range = c(0.01, 0.3),
                                           scale = "linear"))$scale, "linear")
})

test_that("bad bounds are rejected rather than silently accepted", {
  expect_error(.xgb_norm_bound("eta", c(0, 0.3)), "log scale")
  expect_error(.xgb_norm_bound("eta", c(0.3, 0.01)), "below the lower")
  expect_error(.xgb_norm_bound("eta", c(0.1)), "must be a numeric")
})

test_that("eta is sampled log-uniformly, not uniformly", {
  set.seed(7)
  s <- .xgb_sample_bounds(list(eta = c(0.001, 1)), 4000)
  expect_true(all(s$eta >= 0.001 & s$eta <= 1))
  # under a log draw the median sits near the geometric midpoint (~0.0316),
  # whereas a uniform draw would put it near 0.5
  expect_lt(median(s$eta), 0.1)
  # and log(eta) should be close to uniform: mean near the midpoint of the logs
  expect_equal(mean(log(s$eta)), mean(log(c(0.001, 1))), tolerance = 0.1)
})

test_that("integer parameters are integral and alpha has mass at zero", {
  set.seed(7)
  s <- .xgb_sample_bounds(list(max_depth = c(2, 10), alpha = c(0, 5)), 2000)
  expect_true(all(s$max_depth == as.integer(s$max_depth)))
  expect_true(all(s$max_depth >= 2 & s$max_depth <= 10))
  expect_gt(mean(s$alpha == 0), 0.3)   # a substantial share exactly zero
})

test_that("boundary position is reported and flagged", {
  p <- .xgb_bounds_position(list(eta = c(0.01, 0.3)), list(eta = 0.01))
  expect_true(p$at_boundary)
  p2 <- .xgb_bounds_position(list(eta = c(0.01, 0.3)), list(eta = 0.055))
  expect_false(p2$at_boundary)   # geometric middle of a log range
})

## ---- integration: these fit models ---------------------------------------

test_that("the new arguments are inert on the default grid path", {
  ## Seeded deliberately. Without a seed this routine has never been
  ## reproducible: xgboost was left to its own RNG, so two identical calls
  ## disagree. An unseeded comparison here would be a flaky test, not a
  ## compatibility check.
  d <- make_toy()
  args <- list(fixed = y ~ x1 + x2 + x3, smp_data = d, smp_weights = "wt",
               domains = "muni", folds = 3, nround = c(10, 20),
               max_depth = c(2, 3), seed = 101, verbose = FALSE)
  a <- do.call(xgb_tune, args)
  b <- do.call(xgb_tune, c(args, list(search = "grid", n_iter = 999,
                                      bounds = list(eta = c(0.001, 0.9)))))
  expect_equal(a[1:14], b[1:14])
  expect_equal(a$search, "grid")
  expect_equal(a$n_configs, 32L)
  expect_null(a$position)
})

test_that("a seed makes the grid path reproducible, which it was not before", {
  d <- make_toy()
  g <- function() xgb_tune(y ~ x1 + x2 + x3, smp_data = d, smp_weights = "wt",
                           domains = "muni", folds = 3, nround = 10,
                           max_depth = 3, subsample = 0.8, colsample_bytree = 0.6,
                           seed = 202, verbose = FALSE)
  expect_equal(g()[1:14], g()[1:14])
})

test_that("random search draws n_iter configurations and reports position", {
  d <- make_toy()
  r <- xgb_tune(y ~ x1 + x2 + x3, smp_data = d, smp_weights = "wt",
                domains = "muni", folds = 3, search = "random", n_iter = 6,
                bounds = list(eta = c(0.05, 0.3), max_depth = c(2, 4)),
                seed = 99, verbose = FALSE)
  expect_equal(r$search, "random")
  expect_equal(r$n_configs, 6L)
  expect_true(all(c("parameter", "position", "at_boundary") %in% names(r$position)))
  expect_setequal(r$position$parameter, c("eta", "max_depth"))
  # unnamed parameters stay at the first value of their own argument
  expect_equal(r$lambda, 0.5)   # first value of the lambda argument, c(0.5, 1.5)
})

test_that("seed makes the search reproducible", {
  d <- make_toy()
  f <- function() suppressWarnings(xgb_tune(y ~ x1 + x2 + x3, smp_data = d, smp_weights = "wt",
                           domains = "muni", folds = 3, search = "random",
                           n_iter = 4, seed = 123, verbose = FALSE))
  expect_equal(f()[1:14], f()[1:14])
})

test_that("early stopping takes nrounds out of the search", {
  d <- make_toy()
  r <- xgb_tune(y ~ x1 + x2 + x3, smp_data = d, smp_weights = "wt",
                domains = "muni", folds = 3, search = "random", n_iter = 4,
                bounds = list(eta = c(0.05, 0.3)),
                early_stopping_rounds = 5, nrounds_max = 60,
                seed = 5, verbose = FALSE)
  expect_true(!is.null(r$best_iters_by_fold))
  expect_length(r$best_iters_by_fold, 3)
  expect_true(all(r$best_iters_by_fold >= 1))
  # the reported nround is the learned stopping point, not a grid value
  expect_equal(r$nround, round(mean(r$best_iters_by_fold)))
})
