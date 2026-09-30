# megb's two opt-in cross-validation options (6bd056f), both off by default.
# options(megb.cv_folds = "domain"): xgb.cv holds out 5 folds of whole domains, so
#   no domain has rows on both sides of a fold. Default "rows": nfold = 5 random rows.
# options(megb.early_stopping = FALSE): no early stopping; the booster trains to nrounds.
# The tests record the arguments train_gbmodel passes to xgboost::xgb.cv.
.cv_setup <- function() {
  set.seed(11)
  D <- 25; n <- 12
  X <- data.frame(x1 = rnorm(D * n), x2 = runif(D * n))
  g <- rep(sprintf("d%02d", 1:D), each = n)
  y <- 0.5 * X$x1 + rnorm(D, 0, 0.3)[as.integer(factor(g))] + rnorm(D * n)
  list(X = X, y = y, g = g)
}
.record_cv <- function(env) {
  orig <- xgboost::xgb.cv
  function(...) { env$args <- list(...); orig(...) }
}

test_that("default: random-row folds with early stopping, groups ignored", {
  skip_if_not_installed("xgboost")
  S <- .cv_setup(); rec <- new.env()
  local_mocked_bindings(xgb.cv = .record_cv(rec), .package = "xgboost")
  withr::local_options(megb.cv_folds = NULL, megb.early_stopping = NULL)
  r <- povmap:::train_gbmodel("xgboost", S$X, S$y, params = list(eta = 0.3, max_depth = 3, nrounds = 200),
                              groups = S$g)
  expect_equal(rec$args$nfold, 5)
  expect_null(rec$args$folds)
  expect_equal(rec$args$early_stopping_rounds, 10)
  expect_lt(r$best_iter, 200)
})

test_that("megb.cv_folds = 'domain': 5 folds of whole domains covering every row once", {
  skip_if_not_installed("xgboost")
  S <- .cv_setup(); rec <- new.env()
  local_mocked_bindings(xgb.cv = .record_cv(rec), .package = "xgboost")
  withr::local_options(megb.cv_folds = "domain", megb.early_stopping = NULL)
  r <- povmap:::train_gbmodel("xgboost", S$X, S$y, params = list(eta = 0.3, max_depth = 3, nrounds = 200),
                              groups = S$g)
  f <- rec$args$folds
  expect_length(f, 5)
  expect_null(rec$args$nfold)
  expect_setequal(unlist(f), seq_along(S$y))                     # every row held out exactly once
  expect_equal(length(unlist(f)), length(S$y))
  dom_folds <- tapply(rep(seq_along(f), lengths(f)), S$g[unlist(f)], function(k) length(unique(k)))
  expect_true(all(dom_folds == 1L))                               # no domain split across folds
  expect_equal(rec$args$early_stopping_rounds, 10)
  expect_length(r$oof_prediction, length(S$y))
})

test_that("megb.cv_folds = 'domain' without groups falls back to row folds", {
  skip_if_not_installed("xgboost")
  S <- .cv_setup(); rec <- new.env()
  local_mocked_bindings(xgb.cv = .record_cv(rec), .package = "xgboost")
  withr::local_options(megb.cv_folds = "domain")
  povmap:::train_gbmodel("xgboost", S$X, S$y, params = list(eta = 0.3, max_depth = 3, nrounds = 50))
  expect_null(rec$args$folds)
  expect_equal(rec$args$nfold, 5)
})

test_that("megb.early_stopping = FALSE: no early stopping, booster trains to nrounds", {
  skip_if_not_installed("xgboost")
  S <- .cv_setup(); rec <- new.env()
  local_mocked_bindings(xgb.cv = .record_cv(rec), .package = "xgboost")
  withr::local_options(megb.cv_folds = "domain", megb.early_stopping = FALSE)
  r <- povmap:::train_gbmodel("xgboost", S$X, S$y, params = list(eta = 0.3, max_depth = 3, nrounds = 60),
                              groups = S$g)
  expect_null(rec$args$early_stopping_rounds)
  expect_equal(r$best_iter, 60L)
  expect_equal(xgboost::xgb.get.num.boosted.rounds(r$boosting), 60)
})
