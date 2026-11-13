# Convention v0.6 defaults and xgb_tune_repeated().

test_that("the defaults are the convention's", {
  f <- formals(xgb)
  expect_equal(eval(f$area_pool)[1], "out_of_fold"); expect_equal(eval(f$sub_pool)[1], "in_sample")
  expect_false(f$weightedBS); expect_false(f$area_deflate); expect_false(f$sub_deflate)
  expect_equal(formals(xgb_top_configs)$floor, 1); expect_equal(formals(xgb_top_configs)$cap, 8)
  expect_equal(formals(xgb_tune)$n_iter, 150); expect_equal(formals(xgb_tune)$nrounds_max, 3000)
  expect_equal(eval(formals(xgb_tune)$search)[1], "grid")
  b <- .XGB_DEFAULT_BOUNDS
  expect_equal(b$max_depth, c(1, 12)); expect_equal(b$min_child_weight, c(1, 1000)); expect_equal(b$lambda, c(0.01, 10000))
  expect_equal(b$colsample_bytree, c(0.1, 1)); expect_equal(b$nround, c(50, 3000)); expect_equal(b$eta, c(0.01, 0.3)); expect_equal(b$subsample, c(0.5, 1))
})

test_that("xgb_top_configs does not pad the paired one-SE set by default", {
  set.seed(1); n <- 12; k <- 10
  cf <- matrix(rnorm(n * k, sd = 0.01), n, k) + c(0, rep(1, n - 1))      # one clear winner: the set is a single configuration
  grid <- data.frame(nround = 100, max_depth = 3, colsample_bytree = 1, colsample_bylevel = 1, colsample_bynode = 1, subsample = 1,
                     min_child_weight = 1, eta = 0.1, gamma = 0, max_delta_step = 0, lambda = 1, alpha = 0)[rep(1, n), ]
  tn <- list(candidates = grid, cv_by_fold = cf)
  expect_equal(nrow(xgb_top_configs(tn)), 1L)
  expect_equal(nrow(xgb_top_configs(tn, floor = 3)), 3L)
})

test_that("xgb_tune_repeated scores one draw on every fold seed", {
  set.seed(2); d <- 24; per <- 12; muni <- rep(seq_len(d), each = per); x1 <- rnorm(d * per); x2 <- rnorm(d * per)
  dat <- data.frame(y = 1 + x1 - 0.5 * x2 + rnorm(d)[muni] + rnorm(d * per, sd = 0.5), x1 = x1, x2 = x2, muni = muni)
  b <- list(eta = c(0.1, 0.3), max_depth = c(1, 3), nround = c(20, 40))
  r <- suppressWarnings(xgb_tune_repeated(y ~ x1 + x2, dat, domains = "muni", fold_seeds = c(123, 124), n_iter = 5, bounds = b, folds = 3, verbose = FALSE))
  expect_named(r, c("seed_123", "seed_124"))
  expect_equal(r[[1]]$candidates[, c("eta", "max_depth", "nround")], r[[2]]$candidates[, c("eta", "max_depth", "nround")])
  expect_true(all(r[[1]]$candidates$max_depth <= 3 & r[[1]]$candidates$nround >= 20))
  expect_false(isTRUE(all.equal(r[[1]]$cv_by_fold, r[[2]]$cv_by_fold)))
  r2 <- suppressWarnings(xgb_tune_repeated(y ~ x1 + x2, dat, domains = "muni", fold_seeds = c(123, 124), n_iter = 5, bounds = b, folds = 3, verbose = FALSE))
  expect_identical(r[[1]]$cv_by_fold, r2[[1]]$cv_by_fold)
  tc <- xgb_top_configs(r); expect_true(nrow(tc) >= 1 && nrow(tc) <= 8)
})

test_that("a plain xgb() call is the convention: out-of-fold area pool, uniform draws", {
  skip_if_not_installed("xgboost")
  set.seed(7); nd <- 12; ns <- 6
  pop <- expand.grid(sa = seq_len(ns), dom = seq_len(nd)); pop$sub <- sprintf("%02d_%02d", pop$dom, pop$sa); pop$x1 <- rnorm(nrow(pop)); pop$x2 <- runif(nrow(pop))
  pop$npop <- round(runif(nrow(pop), 50, 500)); pop$p <- plogis(-1 + 0.8 * pop$x1 + rnorm(nd, 0, 0.4)[pop$dom]); pop$dom <- sprintf("%02d", pop$dom)
  smp <- do.call(rbind, lapply(sort(sample(unique(pop$dom), 9)), function(d) { s <- pop[pop$dom == d, ][sample(ns, 4), ]; s$y <- rbinom(nrow(s), 25, s$p) / 25; s$wt <- exp(rnorm(nrow(s), 0, 0.5)) * 100; s }))
  fit <- function(...) { r <- NULL
    suppressWarnings(suppressMessages(utils::capture.output(r <- xgb(y ~ x1 + x2, smp, smp_weights = "wt", pop_data = pop, pop_weights = "npop", domains = "dom",
           sub_domains = "sub", transformation = "arcsin", bootstrap = TRUE, B = 4, L = 5, seed = 11, nrounds = 10, max_depth = 2, eta = 0.3, pool_folds = 3, ...))))
    r }
  a <- fit(); b <- fit(area_pool = "out_of_fold", weightedBS = FALSE, sub_pool = "in_sample", area_deflate = FALSE)
  a$out_call <- b$out_call <- NULL; a$model <- b$model <- NULL
  expect_identical(a, b); expect_false(is.null(a$pools))
})
