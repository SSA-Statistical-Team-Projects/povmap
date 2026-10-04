# xgb_tune must fit its candidate models on the weights xgb() fits on.
#
# Until 2026-10 the tuner passed the raw survey weights to xgboost while xgb()
# rescales them within each domain (rescale_weights = TRUE) and divides by their
# mean. With sub-area weights summing to population counts, the tuner's scale
# was over a thousand times the fit's, so min_child_weight almost never bound
# and lambda barely shrank during tuning, then both bit in the fit. In the
# Colombia migration SAE study every tuned xgb configuration was affected.
#
# Tested here, by recording the weight vector xgboost::xgb.train actually
# receives: for each fold, the tuner's fold model gets bit-identical weights to
# xgb() fitted on that fold's training rows, with rescale_weights TRUE and
# FALSE. On the old code the TRUE case fails on the values and the FALSE case
# fails because the argument does not exist.

make_data <- function(seed = 7) {
  set.seed(seed)
  nd <- 30; ns <- 12
  pop <- expand.grid(sa = seq_len(ns), dom = seq_len(nd))
  pop$sub <- sprintf("%02d_%02d", pop$dom, pop$sa)
  pop$x1 <- rnorm(nrow(pop)); pop$x2 <- runif(nrow(pop))
  u <- rnorm(nd, 0, 0.4)[pop$dom]
  pop$npop <- round(runif(nrow(pop), 50, 500))
  pop$p <- plogis(-2 + 0.8 * pop$x1 + 0.6 * pop$x2 + u)
  pop$dom <- sprintf("%02d", pop$dom)
  ## sample: 18 domains, 3 to 9 sub-areas each, survey-scale unequal weights
  sd_ <- sort(sample(unique(pop$dom), 18))
  smp <- do.call(rbind, lapply(sd_, function(d) {
    s <- pop[pop$dom == d, ][sample(ns, sample(3:9, 1)), ]
    s$y <- rbinom(nrow(s), 25, s$p) / 25
    s$wt <- exp(rnorm(nrow(s), 0, 0.7)) * 400
    s }))
  rownames(smp) <- NULL
  list(pop = pop[, c("sub", "dom", "x1", "x2", "npop")],
       smp = smp[, c("sub", "dom", "x1", "x2", "y", "wt")])
}

## every weight vector xgboost::xgb.train receives while expr runs, in call order
capture_train_weights <- function(expr) {
  env <- new.env(); env$w <- list()
  suppressMessages(trace("xgb.train", where = asNamespace("xgboost"), print = FALSE,
    tracer = bquote(assign("w", c(get("w", envir = .(env)),
                                  list(as.numeric(xgboost::getinfo(data, "weight")))),
                           envir = .(env)))))
  on.exit(suppressMessages(untrace("xgb.train", where = asNamespace("xgboost"))))
  force(expr)
  env$w
}

## the tuner's scorer, run in-process (foreach sequential) on fixed folds of domains
tuner_fold_weights <- function(d, fold_of_dom, ...) {
  foreach::registerDoSEQ()
  X_smp <- d$smp[, c("x1", "x2", "dom")]
  cluster_col <- data.frame(dom = d$smp$dom, fold = fold_of_dom[d$smp$dom])
  grid <- data.frame(nround = 3, max_depth = 2, colsample_bytree = 1, colsample_bylevel = 1,
                     colsample_bynode = 1, subsample = 1, min_child_weight = 1, eta = 0.3,
                     gamma = 0, max_delta_step = 0, lambda = 1, alpha = 0)
  capture_train_weights(.xgb_score_configs(
    tunegrid = grid, X_final = X_smp[, c("x1", "x2")], X_smp = X_smp,
    Y_smp = data.frame(labels = d$smp$y), smp_weights = d$smp$wt,
    cluster_col = cluster_col, cluster = "dom", domains = "dom",
    folds = max(fold_of_dom), ...))
}

xgb_fit_weights_on <- function(d, rows, rescale_weights) {
  w <- capture_train_weights(suppressWarnings(povmap::xgb(
    fixed = y ~ x1 + x2, smp_data = d$smp[rows, ], smp_weights = "wt",
    pop_data = d$pop, pop_weights = "npop", domains = "dom", sub_domains = "sub",
    transformation = "no", bootstrap = FALSE, L = 1, seed = 1, nrounds = 3, eta = 0.3,
    max_depth = 2, colsample_bylevel = 1, colsample_bynode = 1,
    rescale_weights = rescale_weights)))
  expect_length(w, 1L)
  w[[1]]
}

fold_of_dom <- function(d, k = 3) {
  doms <- sort(unique(d$smp$dom))
  setNames(rep_len(seq_len(k), length(doms)), doms)
}

test_that("the tuner's fold models get xgb()'s weights (rescale_weights = TRUE, the default)", {
  skip_if_not_installed("xgboost")
  d <- make_data(); f <- fold_of_dom(d)
  wt <- tuner_fold_weights(d, f)                 # default: the tuner must match xgb()'s default
  expect_length(wt, 3L)
  for (k in 1:3) {
    tr <- f[d$smp$dom] != k
    wx <- xgb_fit_weights_on(d, which(tr), rescale_weights = TRUE)
    expect_identical(wt[[k]], wx)
    ## and the scale is the fit's, not the survey's
    expect_equal(mean(wt[[k]]), 1, tolerance = 1e-6)
    expect_gt(mean(d$smp$wt[tr]), 100)
  }
})

test_that("the tuner's fold models get xgb()'s weights (rescale_weights = FALSE)", {
  skip_if_not_installed("xgboost")
  d <- make_data(); f <- fold_of_dom(d)
  wt <- tuner_fold_weights(d, f, rescale_weights = FALSE)
  expect_length(wt, 3L)
  for (k in 1:3) {
    tr <- f[d$smp$dom] != k
    expect_identical(wt[[k]], xgb_fit_weights_on(d, which(tr), rescale_weights = FALSE))
  }
})

test_that("xgb_tune records the weighting it used", {
  skip_if_not_installed("xgboost")
  d <- make_data()
  r <- suppressWarnings(xgb_tune(y ~ x1 + x2, smp_data = d$smp, smp_weights = "wt",
                                 domains = "dom", folds = 3, nround = 3, max_depth = 2,
                                 colsample_bytree = 1, min_child_weight = 1, lambda = 1,
                                 seed = 1, cpus = 1, verbose = FALSE))
  expect_true(isTRUE(r$rescale_weights))
  r2 <- suppressWarnings(xgb_tune(y ~ x1 + x2, smp_data = d$smp, smp_weights = "wt",
                                  domains = "dom", folds = 3, nround = 3, max_depth = 2,
                                  colsample_bytree = 1, min_child_weight = 1, lambda = 1,
                                  seed = 1, cpus = 1, verbose = FALSE, rescale_weights = FALSE))
  expect_false(r2$rescale_weights)
  ## the fold assignment uses no weights, so it is the same either way
  expect_identical(r$fold_assignment, r2$fold_assignment)
})

## variance_y: xgb() multiplies its fit's weights by variance_y^-0.5 and divides by
## their mean once more (point_estim_xgb); the tuner must pass the same weights.
with_variance <- function(d, seed = 11) {
  set.seed(seed)
  d$smp$v <- runif(nrow(d$smp), 0.2, 3)
  d
}

xgb_fit_weights_var_on <- function(d, rows, rescale_weights) {
  w <- capture_train_weights(suppressWarnings(povmap::xgb(
    fixed = y ~ x1 + x2, smp_data = d$smp[rows, ], smp_weights = "wt",
    pop_data = d$pop, pop_weights = "npop", domains = "dom", sub_domains = "sub",
    transformation = "no", bootstrap = FALSE, L = 1, seed = 1, nrounds = 3, eta = 0.3,
    max_depth = 2, colsample_bylevel = 1, colsample_bynode = 1,
    rescale_weights = rescale_weights, variance_y = "v")))
  expect_length(w, 1L)
  w[[1]]
}

test_that("with variance_y, the tuner's fold models get xgb()'s weights", {
  skip_if_not_installed("xgboost")
  d <- with_variance(make_data()); f <- fold_of_dom(d)
  for (rw in c(TRUE, FALSE)) {
    ## the factor exactly as xgb_tune(variance_y = "v") computes it
    wt <- tuner_fold_weights(d, f, rescale_weights = rw, het = as.numeric(d$smp[["v"]])^-0.5)
    expect_length(wt, 3L)
    for (k in 1:3) {
      tr <- f[d$smp$dom] != k
      wx <- xgb_fit_weights_var_on(d, which(tr), rescale_weights = rw)
      expect_identical(wt[[k]], wx, info = paste("rescale_weights =", rw, "fold", k))
      ## the factor matters: these are not the weights without variance_y
      expect_false(isTRUE(all.equal(wx, xgb_fit_weights_on(d, which(tr), rescale_weights = rw))))
    }
  }
})

test_that("xgb_tune(variance_y =) scores with those weights and records variance_y", {
  skip_if_not_installed("xgboost")
  d <- with_variance(make_data())
  args <- list(y ~ x1 + x2, smp_data = d$smp, smp_weights = "wt", domains = "dom", folds = 3,
               nround = 3, max_depth = 2, colsample_bytree = 1, subsample = 1,
               min_child_weight = 1, lambda = 1, seed = 1, cpus = 1, verbose = FALSE)
  r  <- suppressWarnings(do.call(xgb_tune, c(args, list(variance_y = "v"))))
  r0 <- suppressWarnings(do.call(xgb_tune, args))
  expect_identical(r$variance_y, "v")
  expect_null(r0$variance_y)
  expect_identical(r$fold_assignment, r0$fold_assignment)
  expect_false(isTRUE(all.equal(r$cv_by_fold, r0$cv_by_fold)))
  ## the same scores as the scorer given variance_y^-0.5 directly, on the same folds
  foreach::registerDoSEQ()
  fa <- r$fold_assignment
  X_smp <- d$smp[, c("x1", "x2", "dom")]
  sc <- .xgb_score_configs(
    tunegrid = r$candidates[1, names(r$candidates) %in% c("nround", "max_depth", "colsample_bytree",
      "colsample_bylevel", "colsample_bynode", "subsample", "min_child_weight", "eta", "gamma",
      "max_delta_step", "lambda", "alpha")],
    X_final = X_smp[, c("x1", "x2")], X_smp = X_smp, Y_smp = data.frame(labels = d$smp$y),
    smp_weights = d$smp$wt, cluster_col = data.frame(dom = d$smp$dom, fold = fa$fold[match(d$smp$dom, fa$dom)]),
    cluster = "dom", domains = "dom", folds = 3, seed = 1, het = as.numeric(d$smp$v)^-0.5)
  expect_equal(unname(sc$OPT), unname(r$cv_by_fold), tolerance = 1e-12)
})

test_that("variance_y must name a positive numeric variable", {
  skip_if_not_installed("xgboost")
  d <- with_variance(make_data())
  run <- function(smp, vy) suppressWarnings(xgb_tune(y ~ x1 + x2, smp_data = smp, smp_weights = "wt",
    domains = "dom", folds = 3, nround = 3, max_depth = 2, colsample_bytree = 1,
    min_child_weight = 1, lambda = 1, seed = 1, cpus = 1, verbose = FALSE, variance_y = vy))
  expect_error(run(d$smp, "nope"), "variance_y must name")
  s <- d$smp; s$v[1] <- 0
  expect_error(run(s, "v"), "positive")
  s <- d$smp; s$v[2] <- NA
  expect_error(run(s, "v"), "positive")
})
