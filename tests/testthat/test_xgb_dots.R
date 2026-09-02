# `...` in xgb_tune() is forwarded to xgboost::xgb.train(). Anything that is
# neither a formal of xgb_tune() nor an xgboost parameter is therefore
# discarded silently. That is how a caller written against a newer signature
# runs against an older library and returns a plausible-looking result: when
# `search`, `n_iter`, `bounds` and `seed` were added in 3b0cb07, supplying them
# to the build that preceded it evaluated the 32-point default grid -- every
# point at eta = 0.3 -- and returned it as if it were a random search.
#
# These tests pin the guard that makes that loud.

test_that("an unrecognised name in ... is an error, not a silent discard", {
  expect_error(.check_xgb_dots(list(search = "random")), "unrecognised argument")
  expect_error(.check_xgb_dots(list(n_iter = 150)),      "unrecognised argument")
  expect_error(.check_xgb_dots(list(bounds = list())),   "unrecognised argument")
})

test_that("the message names the arguments and explains the consequence", {
  e <- tryCatch(.check_xgb_dots(list(search = "random", n_iter = 150)),
                error = function(e) conditionMessage(e))
  expect_match(e, "'search'")
  expect_match(e, "'n_iter'")
  expect_match(e, "silently discarded")
  expect_match(e, "NEWER xgb_tune")
})

test_that("a near-miss spelling of a real formal is named", {
  e <- tryCatch(.check_xgb_dots(list(fold = 5), fn_formals = c("folds", "nround", "eta")),
                error = function(e) conditionMessage(e))
  expect_match(e, "did you mean")
  expect_match(e, "folds")
})

test_that("genuine xgboost parameters still pass through", {
  expect_true(.check_xgb_dots(list(objective = "reg:squarederror")))
  expect_true(.check_xgb_dots(list(tree_method = "hist", max_bin = 256)))
  expect_true(.check_xgb_dots(list(nthread = 4, eval_metric = "rmse")))
  expect_true(.check_xgb_dots(list()))
})

test_that("an unnamed argument in ... is an error", {
  expect_error(.check_xgb_dots(list(5)), "without a name")
  expect_error(.check_xgb_dots(stats::setNames(list(1, 2), c("objective", ""))), "without a name")
})

test_that("the escape hatch exists so a new upstream parameter cannot be blocked", {
  expect_true(.check_xgb_dots(list(brand_new_xgb_param = 1),
                              allow = "brand_new_xgb_param"))
  expect_warning(.check_xgb_dots(list(brand_new_xgb_param = 1), strict = FALSE),
                 "unrecognised argument")
})

test_that("xgb_tune() itself rejects an argument its signature does not have", {
  skip_if_not_installed("xgboost")
  set.seed(1)
  d <- data.frame(y = rnorm(200), x1 = rnorm(200), x2 = rnorm(200),
                  muni = rep(1:20, each = 10), wt = 1)
  # `not_an_argument` must fail BEFORE any model is fitted, so this is fast
  expect_error(
    xgb_tune(fixed = y ~ x1 + x2, smp_data = d, smp_weights = "wt",
             domains = "muni", folds = 2, cpus = 1, verbose = FALSE,
             not_an_argument = TRUE),
    "unrecognised argument")
})
