# megb's internal cross-validation settings: arguments cv_folds and early_stopping
# of megb() (and of train_gbmodel, where they act), with the global options
# megb.cv_folds / megb.early_stopping as the fallback when an argument is omitted.
# cv_folds = "domain" (default): xgb.cv holds out 5 folds of whole domains, so no
#   domain has rows on both sides of a fold. "rows": nfold = 5 random rows.
#   "domain" with no domain labels falls back to "rows" with a message.
# early_stopping = TRUE (default): stop after 10 rounds without improvement;
#   FALSE: fold models and booster train to nrounds.
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
.fit <- function(S, groups = S$g, nrounds = 200, ...)
  povmap:::train_gbmodel("xgboost", S$X, S$y, params = list(eta = 0.3, max_depth = 3, nrounds = nrounds),
                         groups = groups, ...)
## K: the number of folds expected, 10 by default from povmap 1.0.1 (it was 5)
.is_domain_folds <- function(args, g, K = 10L) {
  f <- args$folds
  if (is.null(f) || length(f) != K || !is.null(args$nfold)) return(FALSE)
  if (!setequal(unlist(f), seq_along(g)) || length(unlist(f)) != length(g)) return(FALSE)
  all(tapply(rep(seq_along(f), lengths(f)), g[unlist(f)], function(k) length(unique(k))) == 1L)
}
.is_row_folds <- function(args, K = 10) is.null(args$folds) && identical(args$nfold, K)

test_that("default, no argument and no option: domain folds with early stopping", {
  skip_if_not_installed("xgboost")
  S <- .cv_setup(); rec <- new.env()
  local_mocked_bindings(xgb.cv = .record_cv(rec), .package = "xgboost")
  withr::local_options(megb.cv_folds = NULL, megb.early_stopping = NULL)
  r <- .fit(S)
  expect_true(.is_domain_folds(rec$args, S$g))
  expect_equal(rec$args$early_stopping_rounds, 10)
  expect_lt(r$best_iter, 200)
  expect_length(r$oof_prediction, length(S$y))
})

test_that("cv_folds = 'rows': folds of random rows (ten by default)", {
  skip_if_not_installed("xgboost")
  S <- .cv_setup(); rec <- new.env()
  local_mocked_bindings(xgb.cv = .record_cv(rec), .package = "xgboost")
  withr::local_options(megb.cv_folds = NULL, megb.early_stopping = NULL)
  .fit(S, cv_folds = "rows")
  expect_true(.is_row_folds(rec$args))
  expect_equal(rec$args$early_stopping_rounds, 10)
})

test_that("with the argument omitted, the global options still apply", {
  skip_if_not_installed("xgboost")
  S <- .cv_setup(); rec <- new.env()
  local_mocked_bindings(xgb.cv = .record_cv(rec), .package = "xgboost")
  withr::local_options(megb.cv_folds = "rows", megb.early_stopping = FALSE)
  r <- .fit(S, nrounds = 40)
  expect_true(.is_row_folds(rec$args))
  expect_null(rec$args$early_stopping_rounds)
  expect_equal(r$best_iter, 40L)
})

test_that("an explicit argument overrides the global option", {
  skip_if_not_installed("xgboost")
  S <- .cv_setup(); rec <- new.env()
  local_mocked_bindings(xgb.cv = .record_cv(rec), .package = "xgboost")
  withr::local_options(megb.cv_folds = "rows", megb.early_stopping = FALSE)
  .fit(S, cv_folds = "domain", early_stopping = TRUE)
  expect_true(.is_domain_folds(rec$args, S$g))
  expect_equal(rec$args$early_stopping_rounds, 10)
  withr::local_options(megb.cv_folds = "domain", megb.early_stopping = TRUE)
  .fit(S, cv_folds = "rows", early_stopping = FALSE, nrounds = 40)
  expect_true(.is_row_folds(rec$args))
  expect_null(rec$args$early_stopping_rounds)
})

test_that("cv_folds = 'domain' without domain labels falls back to row folds with a message", {
  skip_if_not_installed("xgboost")
  S <- .cv_setup(); rec <- new.env()
  local_mocked_bindings(xgb.cv = .record_cv(rec), .package = "xgboost")
  withr::local_options(megb.cv_folds = NULL)
  expect_message(.fit(S, groups = NULL, nrounds = 50), "no domain labels")
  expect_true(.is_row_folds(rec$args))
  expect_silent(suppressWarnings(.fit(S, groups = NULL, nrounds = 50, cv_folds = "rows")))
})

test_that("early_stopping = FALSE: no early stopping, booster trains to nrounds", {
  skip_if_not_installed("xgboost")
  S <- .cv_setup(); rec <- new.env()
  local_mocked_bindings(xgb.cv = .record_cv(rec), .package = "xgboost")
  withr::local_options(megb.cv_folds = NULL, megb.early_stopping = NULL)
  r <- .fit(S, nrounds = 60, early_stopping = FALSE)
  expect_null(rec$args$early_stopping_rounds)
  expect_equal(r$best_iter, 60L)
  expect_equal(xgboost::xgb.get.num.boosted.rounds(r$boosting), 60)
})

# The same through megb(): its default is domain folds, and its arguments reach
# train_gbmodel and override the options.
.megb_data <- function() {
  set.seed(3)
  D <- 30; n_pop <- 20; n_smp <- 8
  pop <- data.frame(dom = factor(rep(sprintf("d%02d", 1:D), each = n_pop)),
                    x1 = rnorm(D * n_pop), x2 = runif(D * n_pop), w = rpois(D * n_pop, 20) + 1)
  pop$y <- 0.5 * pop$x1 + rnorm(D, 0, 0.3)[as.integer(pop$dom)] + rnorm(nrow(pop))
  smp <- do.call(rbind, lapply(split(pop, pop$dom)[1:25], function(d) d[sample(nrow(d), n_smp), ]))
  smp$wt <- 1
  list(pop = pop, smp = smp)
}
.megb_fit <- function(M, ...) suppressWarnings(povmap::megb(
  fixed = y ~ x1 + x2, smp_data = M$smp, smp_weights = "wt", pop_data = M$pop[, c("dom", "x1", "x2", "w")],
  pop_weights = "w", domains = "dom", transformation = "no", mse = FALSE, na.rm = FALSE, seed = 5,
  gradient_params = list(eta = 0.3, max_depth = 3, subsample = 1, nrounds = 100), ...))

test_that("megb(): default is domain folds; its arguments override the options", {
  skip_if_not_installed("xgboost")
  M <- .megb_data(); rec <- new.env()
  local_mocked_bindings(xgb.cv = .record_cv(rec), .package = "xgboost")
  orig_tg <- povmap:::train_gbmodel                 # record the domain labels megb passes, in its row order
  local_mocked_bindings(train_gbmodel = function(...) { rec$groups <- list(...)$groups; orig_tg(...) },
                        .package = "povmap")
  withr::local_options(megb.cv_folds = NULL, megb.early_stopping = NULL)
  suppressMessages(.megb_fit(M))
  expect_setequal(as.character(rec$groups), as.character(unique(M$smp$dom)))
  expect_true(.is_domain_folds(rec$args, as.character(rec$groups)))
  expect_equal(rec$args$early_stopping_rounds, 10)
  withr::local_options(megb.cv_folds = "domain", megb.early_stopping = TRUE)
  suppressMessages(.megb_fit(M, cv_folds = "rows", early_stopping = FALSE))
  expect_true(.is_row_folds(rec$args))
  expect_null(rec$args$early_stopping_rounds)
  expect_error(.megb_fit(M, cv_folds = "municipality"))
  expect_error(.megb_fit(M, early_stopping = NA))
})
