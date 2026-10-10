# Least-overfit configuration selection: xgb_top_configs(rule = "least_overfit" / "one_se_ratio"),
# xgb_overfit_ratio() and xgb_tune(overfit_ratio = ).
#
# The rule (DRC, 2026-10-08): among the tuning candidates in the paired one-SE set, take the single
# configuration with the lowest out-of-fold / in-sample residual ratio. The DRC fixture test checks
# the selection step against the 16 indicators of the DRC application. Everything is opt-in: the
# default outputs of xgb_top_configs() and xgb_tune() are checked unchanged.

## A synthetic tune whose paired one-SE set is known by construction. Candidate i has fold scores
## m[i] + a[i] * e with e = (-1, 1, -1, 1); the best candidate (lowest m) has a = 0. Its fold-wise
## difference from the best is then d_i = m_i - m_best on average with standard error
## a_i * sd(e) / 2 = a_i * 0.57735, so i is in the set exactly when a_i * 0.57735 >= d_i.
synth_tune <- function(m, a, extra = list()) {
  e <- c(-1, 1, -1, 1)
  cv <- outer(m, rep(1, 4)) + outer(a, e)
  n <- length(m)
  cand <- data.frame(nround = 100 + seq_len(n), max_depth = 3, colsample_bytree = 0.8, colsample_bylevel = 1,
                     colsample_bynode = 1, subsample = 0.8, min_child_weight = 5, eta = 0.3, gamma = 0,
                     max_delta_step = 0, lambda = 1, alpha = 0)
  cand$mse <- rowMeans(cv); cand$rank <- rank(cand$mse)
  c(list(candidates = cand, cv_by_fold = cv), extra)
}
m5 <- c(1.0, 1.1, 1.2, 1.3, 1.4)
a5 <- c(0, 0.2, 0.4, 0.1, 0.1)             # in the set: rows 1, 2, 3 (se 0.115, 0.231 vs d 0.1, 0.2); 4 and 5 out

test_that("synthetic tune: the one-SE set is rows 1 to 3 (rule is the existing paired rule)", {
  t5 <- synth_tune(m5, a5)
  p <- xgb_top_configs(t5, floor = 1, cap = 8)
  expect_equal(p$nrounds, 101:103)
  expect_false("overfit_ratio" %in% names(p))
})

test_that("least_overfit picks the lowest ratio inside the one-SE set only", {
  t5 <- synth_tune(m5, a5)
  r <- c(2.0, 1.5, 1.7, 0.5, 0.4)            # rows 4 and 5 have the lowest ratios but are outside the set
  p <- xgb_top_configs(t5, rule = "least_overfit", ratio = r)
  expect_equal(nrow(p), 1L)
  expect_equal(p$nrounds, 102)
  expect_equal(p$overfit_ratio, 1.5)
  expect_equal(p$weight, 1)
  expect_equal(p$rank, 2L)
  expect_equal(p$cv_score, mean(m5[2] + a5[2] * c(-1, 1, -1, 1)))
  ## floor and cap do not apply to the single-configuration rule
  expect_identical(xgb_top_configs(t5, rule = "least_overfit", ratio = r, floor = 5, cap = 5), p)
})

test_that("least_overfit: ties go to the better score; non-finite ratios cannot be chosen", {
  t5 <- synth_tune(m5, a5)
  expect_equal(xgb_top_configs(t5, rule = "least_overfit", ratio = c(1.5, 1.5, 1.5, 0, 0))$nrounds, 101)
  expect_equal(xgb_top_configs(t5, rule = "least_overfit", ratio = c(1.5, 1.4, 1.4, 0, 0))$nrounds, 102)
  expect_equal(xgb_top_configs(t5, rule = "least_overfit", ratio = c(NA, 1.4, 1.2, 0, 0))$nrounds, 103)
  expect_equal(xgb_top_configs(t5, rule = "least_overfit", ratio = c(NA, Inf, 1.2, 0, 0))$nrounds, 103)
  expect_error(xgb_top_configs(t5, rule = "least_overfit", ratio = c(NA, NA, NA, 1, 1)), "finite ratio")
  ## candidates tied on score: grid order decides which is better, and then which wins a ratio tie
  tt <- synth_tune(c(1, 1, 1.5), c(0, 0, 0))
  tt$cv_by_fold[2, ] <- tt$cv_by_fold[1, ] + c(-.1, .1, -.1, .1)   # same mean, noisier: in the set
  expect_equal(xgb_top_configs(tt, rule = "least_overfit", ratio = c(1, 1, 1))$nrounds, 101)
})

test_that("ratio_cap limits the members considered to the best by score", {
  t5 <- synth_tune(m5, a5)
  r <- c(2.0, 1.5, 1.7, 0.5, 0.4)
  expect_equal(xgb_top_configs(t5, rule = "least_overfit", ratio = r, ratio_cap = 1)$nrounds, 101)
  expect_equal(xgb_top_configs(t5, rule = "least_overfit", ratio = r, ratio_cap = 2)$nrounds, 102)
  expect_equal(xgb_top_configs(t5, rule = "least_overfit", ratio = r, ratio_cap = NULL)$nrounds, 102)
  ## a member beyond the cap with the lowest ratio is not chosen
  expect_equal(xgb_top_configs(t5, rule = "least_overfit", ratio = c(2, 1.5, 0.1, 9, 9), ratio_cap = 2)$nrounds, 102)
  expect_equal(xgb_top_configs(t5, rule = "least_overfit", ratio = c(2, 1.5, 0.1, 9, 9), ratio_cap = 3)$nrounds, 103)
})

test_that("one_se_ratio averages the set members at or below max_ratio, with the fallback", {
  t5 <- synth_tune(m5, a5)
  r <- c(2.0, 1.5, 1.7, 0.5, 0.4)
  p <- xgb_top_configs(t5, rule = "one_se_ratio", ratio = r, max_ratio = 1.7)      # <= includes the boundary
  expect_equal(p$nrounds, c(102, 103))
  expect_equal(p$weight, c(0.5, 0.5))
  expect_equal(p$overfit_ratio, c(1.5, 1.7))
  expect_equal(xgb_top_configs(t5, rule = "one_se_ratio", ratio = r, max_ratio = 5)$nrounds, 101:103)  # whole set
  ## cap limits the number averaged, best scores first
  expect_equal(xgb_top_configs(t5, rule = "one_se_ratio", ratio = r, max_ratio = 5, cap = 2)$nrounds, 101:102)
  ## none qualifies: the single lowest-ratio member (not NA, not outside the set)
  p0 <- xgb_top_configs(t5, rule = "one_se_ratio", ratio = r, max_ratio = 1.0)
  expect_equal(p0$nrounds, 102); expect_equal(p0$weight, 1)
  ## exactly one qualifies: that one (it is the lowest)
  p1 <- xgb_top_configs(t5, rule = "one_se_ratio", ratio = c(1.4, 1.6, 2.5, 0.1, 3), max_ratio = 1.5)
  expect_equal(p1$nrounds, 101); expect_equal(nrow(p1), 1L)
  expect_equal(xgb_top_configs(t5, rule = "one_se_ratio", ratio = c(1.4, 1.6, 2.5, 0.1, 3), max_ratio = 1.7)$nrounds, 101:102)
  ## the default threshold is 2
  expect_equal(xgb_top_configs(t5, rule = "one_se_ratio", ratio = c(1.9, 2.0, 2.1, 0, 0))$nrounds, 101:102)
  ## the table is accepted by xgb(configs = ) as it is
  chk <- povmap:::.xgb_check_configs(p, defaults = povmap:::.xgb_xgb_defaults())
  expect_equal(chk$nrounds, c(102, 103)); expect_equal(chk$weight, c(0.5, 0.5))
})

test_that("ratios are read from the tune objects and averaged over those that have them", {
  r1 <- c(2.0, 1.5, 1.7, 0.5, 0.4); r2 <- c(2.0, 1.9, 1.1, 0.5, 0.4)
  t1 <- synth_tune(m5, a5, list(overfit = data.frame(ratio_total = r1)))
  t2 <- synth_tune(m5, a5, list(overfit = data.frame(ratio_total = r2)))
  t3 <- synth_tune(m5, a5)
  expect_equal(xgb_top_configs(t1, rule = "least_overfit")$nrounds, 102)
  ## averaged over t1 and t2 (t3 has none): row 2 (1.5 + 1.9) / 2 = 1.7, row 3 (1.7 + 1.1) / 2 = 1.4
  pm <- xgb_top_configs(list(t1, t2, t3), rule = "least_overfit")
  expect_equal(pm$nrounds, 103); expect_equal(pm$overfit_ratio, 1.4)
  expect_error(xgb_top_configs(t3, rule = "least_overfit"), "overfit_ratio = TRUE")
  expect_error(xgb_top_configs(t1, rule = "least_overfit", ratio = 1:3), "one element per candidate")
})

test_that("defaults are unchanged: no new column, same set, same rule as before", {
  t5 <- synth_tune(m5, a5)
  p <- xgb_top_configs(t5)                               # floor 3, cap 8, one_se_paired
  expect_equal(p$nrounds, 101:103)
  expect_identical(names(p), c("nrounds", "max_depth", "colsample_bytree", "colsample_bylevel", "colsample_bynode",
                               "subsample", "min_child_weight", "eta", "gamma", "max_delta_step", "lambda", "alpha",
                               "weight", "cv_score", "rank"))
  expect_identical(p, xgb_top_configs(t5, rule = "one_se_paired"))
  ## a ratio on the tune does not alter the other rules
  t5r <- synth_tune(m5, a5, list(overfit = data.frame(ratio_total = c(2, 1.5, 1.7, .5, .4))))
  expect_identical(xgb_top_configs(t5r), p)
  expect_equal(xgb_top_configs(t5, rule = "k_best", k = 2, floor = 1)$nrounds, 101:102)
})

## ---- xgb_overfit_ratio and xgb_tune(overfit_ratio = TRUE) ------------------------------------------

make_ratio_data <- function(seed = 3) {
  set.seed(seed)
  nd <- 30; per <- 5
  d <- data.frame(dom = rep(sprintf("d%02d", seq_len(nd)), each = per))
  d$x1 <- rnorm(nrow(d)); d$x2 <- runif(nrow(d))
  d$y <- plogis(-1 + 0.8 * d$x1 + 0.6 * d$x2 + rnorm(nd, 0, 0.4)[as.integer(factor(d$dom))] + rnorm(nrow(d), 0, 0.5))
  d$wt <- exp(rnorm(nrow(d), 0, 0.5)) * 100
  d$v <- runif(nrow(d), 0.5, 2)
  d
}
ratio_grid <- data.frame(nround = c(15, 30, 20), max_depth = c(2, 3, 2), eta = c(0.3, 0.3, 0.1),
                         min_child_weight = c(1, 1, 3), lambda = c(1, 1, 5), colsample_bytree = 0.8,
                         colsample_bylevel = 1, colsample_bynode = 1, subsample = 0.8, gamma = 0,
                         max_delta_step = 0, alpha = 0)

test_that("xgb_overfit_ratio matches a hand computation with plain xgboost", {
  skip_if_not_installed("xgboost")
  d <- make_ratio_data()
  rr <- suppressWarnings(xgb_overfit_ratio(y ~ x1 + x2, d, ratio_grid, smp_weights = "wt", domains = "dom",
                                           transformation = "arcsin", folds = 4, seed = 901, variance_y = "v"))
  expect_equal(rr$row, 1:3)
  ## by hand, row 2
  yt <- asin(sqrt(d$y)); X <- data.matrix(d[, c("x1", "x2")]); het <- d$v^-0.5
  set.seed(901); cu <- unique(d$dom); fold <- sample(1:4, length(cu), replace = TRUE)[match(d$dom, cu)]
  g <- ratio_grid[2, ]
  p <- list(max_depth = g$max_depth, colsample_bytree = g$colsample_bytree, colsample_bylevel = 1, colsample_bynode = 1,
            subsample = g$subsample, min_child_weight = g$min_child_weight, eta = g$eta, gamma = 0, max_delta_step = 0,
            lambda = g$lambda, alpha = 0, seed = 123, nthread = 1)
  fw <- function(tr) { w <- d$wt[tr]; w <- w * ave(w, d$dom[tr], FUN = length) / ave(w, d$dom[tr], FUN = sum)
    w <- w / mean(w); w <- w * het[tr] / mean(w * het[tr]); w }
  fitm <- function(tr) xgboost::xgb.train(xgboost::xgb.DMatrix(X[tr, ], label = yt[tr], weight = fw(tr)), params = p,
                                          nrounds = g$nround, verbose = 0)
  sdw <- function(x, w) sqrt(sum(w * (x - sum(w * x) / sum(w))^2) / sum(w))
  rin <- yt - predict(fitm(rep(TRUE, nrow(d))), X)
  pr <- numeric(nrow(d)); for (k in 1:4) pr[fold == k] <- predict(fitm(fold != k), X[fold == k, ])
  expect_equal(rr$in_total[2], sdw(rin, d$wt), tolerance = 1e-10)
  expect_equal(rr$oof_total[2], sdw(yt - pr, d$wt), tolerance = 1e-10)
  expect_equal(rr$ratio_total[2], sdw(yt - pr, d$wt) / sdw(rin, d$wt), tolerance = 1e-10)
  expect_true(all(rr$ratio_total > 1))                       # held-out clusters are never fitted better than the sample
  ## rows subset, and nrounds accepted for nround
  r2 <- suppressWarnings(xgb_overfit_ratio(y ~ x1 + x2, d, setNames(ratio_grid, sub("^nround$", "nrounds", names(ratio_grid))),
                                           smp_weights = "wt", domains = "dom", transformation = "arcsin", folds = 4,
                                           seed = 901, variance_y = "v", rows = 2L))
  expect_equal(r2$row, 2L); expect_equal(r2$ratio_total, rr$ratio_total[2])
})

test_that("xgb_tune(overfit_ratio = TRUE) adds the ratios and changes nothing else", {
  skip_if_not_installed("xgboost")
  skip_on_cran()
  d <- make_ratio_data()
  run <- function(...) suppressWarnings(suppressMessages(xgb_tune(
    y ~ x1 + x2, smp_data = d, smp_weights = "wt", domains = "dom", transformation = "arcsin", folds = 4,
    grid = ratio_grid, seed = 5, variance_y = "v", verbose = FALSE, ...)))
  t0 <- run(); t1 <- run(overfit_ratio = TRUE)
  expect_null(t0$overfit)
  t1b <- t1; t1b$overfit <- NULL
  expect_identical(t1b, t0)
  expect_equal(nrow(t1$overfit), 3L)
  ## the same ratios as xgb_overfit_ratio() on the same folds (seed 5) and fit seed
  rr <- suppressWarnings(xgb_overfit_ratio(y ~ x1 + x2, d, ratio_grid, smp_weights = "wt", domains = "dom",
                                           transformation = "arcsin", folds = 4, seed = 5, fit_seed = 5,
                                           variance_y = "v"))
  expect_equal(t1$overfit, rr, tolerance = 1e-10)
  ## and xgb_top_configs reads them
  pk <- xgb_top_configs(t1, rule = "least_overfit", ratio_cap = NULL)
  expect_equal(nrow(pk), 1L)
  expect_true(pk$overfit_ratio %in% t1$overfit$ratio_total)
})

## ---- the DRC selection: 16 indicators --------------------------------------------------------------
## The fixture holds, for each indicator, the candidates the DRC selection computed ratios for (members of the
## one-SE set among the 40 best by fresh score, plus the earlier selection and the best-scoring one), with
## their fresh-fold score, total ratio and one-SE-set flag, and the configuration the DRC selected
## (xgb_all_tunes_least_overfit). The per-fold fresh scores (seeds 901, 902) were not saved, so the set
## itself cannot be recomputed from the fixture: each indicator's tune object is built so that its paired
## one-SE set is exactly the DRC set and its scores are the DRC's fresh scores, and the selection step is
## run through xgb_top_configs().
drc_fix_path <- testthat::test_path("fixtures", "least_overfit_drc_fixture.rds")

drc_synth_tune <- function(x) {
  cd <- x$cand; cd <- cd[order(cd$fresh_mse, cd$row), ]       # best first, ties by grid row as in the DRC
  e <- c(-1, 1, -1, 1)
  s_e <- sqrt(stats::var(e) / 4)
  best <- cd[1, ]
  stopifnot(best$in_set)
  a <- ifelse(cd$in_set, 1.5 * (cd$fresh_mse - best$fresh_mse) / s_e + 1e-6, 0)
  a[1] <- 0
  cv <- outer(cd$fresh_mse, rep(1, 4)) + outer(a, e)
  cand <- cd[, c("nround", "max_depth", "colsample_bytree", "colsample_bylevel", "colsample_bynode", "subsample",
                 "min_child_weight", "eta", "gamma", "max_delta_step", "lambda", "alpha")]
  cand$mse <- rowMeans(cv)
  list(candidates = cand, cv_by_fold = cv, ratio = cd$ratio_total)
}

test_that("DRC set B: least_overfit reproduces the saved selection for all 16 indicators", {
  skip_if(!file.exists(drc_fix_path), "DRC fixture not available")
  fix <- readRDS(drc_fix_path)
  expect_length(fix, 16)
  hp <- c("nround", "max_depth", "colsample_bytree", "colsample_bylevel", "colsample_bynode", "subsample",
          "min_child_weight", "eta", "gamma", "max_delta_step", "lambda", "alpha")
  for (ind in names(fix)) {
    x <- fix[[ind]]
    tn <- drc_synth_tune(x)
    ## the constructed one-SE set is the DRC set: the paired rule returns exactly the flagged rows
    ps <- xgb_top_configs(tn[c("candidates", "cv_by_fold")], floor = 1, cap = nrow(tn$candidates))
    expect_equal(nrow(ps), sum(x$cand$in_set), info = ind)
    p <- xgb_top_configs(tn[c("candidates", "cv_by_fold")], rule = "least_overfit", ratio = tn$ratio)
    expect_equal(nrow(p), 1L, info = ind)
    got <- p[1, ]; names(got)[names(got) == "nrounds"] <- "nround"
    expect_identical(as.numeric(got[1, hp]), as.numeric(x$expected[1, hp]), info = ind)
    ## it is the lowest-ratio member (first in score order on ties), and was ratioed
    mem <- x$cand[x$cand$in_set, ]
    expect_equal(p$overfit_ratio, min(mem$ratio_total), info = ind)
  }
})

test_that("DRC fixture agrees with the saved DRC outputs when they are available", {
  base <- "H:/david/SAE/DRC/data/candidate_estimates/least_overfit_20261008"
  skip_on_cran()
  skip_if(!file.exists(file.path(base, "xgb_all_tunes_least_overfit")) || !file.exists(drc_fix_path),
          "DRC saved outputs not available")
  fix <- readRDS(drc_fix_path)
  tun <- readRDS(file.path(base, "xgb_all_tunes_least_overfit"))
  sel <- utils::read.csv(file.path(base, "least_overfit_selection.csv"))
  mem <- utils::read.csv(file.path(base, "least_overfit_all_members.csv"))
  hp <- c("nround", "max_depth", "colsample_bytree", "colsample_bylevel", "colsample_bynode", "subsample",
          "min_child_weight", "eta", "gamma", "max_delta_step", "lambda", "alpha")
  expect_setequal(names(fix), names(tun))
  for (ind in names(fix)) {
    expect_identical(as.numeric(unlist(tun[[ind]][hp])), as.numeric(fix[[ind]]$expected[1, hp]), info = ind)
    m <- mem[mem$indicator == ind, ]
    expect_equal(fix[[ind]]$cand$ratio_total, m$ratio_total, info = ind)
    expect_equal(sum(fix[[ind]]$cand$in_set), sum(m$in_set), info = ind)
    expect_equal(min(m$ratio_total[m$in_set]), sel$sel_ratio_total[sel$indicator == ind], tolerance = 1e-12, info = ind)
  }
})
