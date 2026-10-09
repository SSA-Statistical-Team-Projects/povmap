# Tests for xgb(interval_scale = "transformed") and xgb(keep_replicates = TRUE).
#
# The natural-scale benchmarked interval is the estimate plus quantiles of the centred bootstrap error
# (benchmarked replicate minus simulated truth) on the scale of the estimates, truncated to [0, 1] for
# arcsin. Near a bound the error distribution is cut off, so many intervals end exactly on 0 or 1, and a
# log outcome (an inequality index) gets negative lower bounds because nothing truncates it.
# interval_scale = "transformed" takes the same quantiles on the transformation scale and back-transforms.
#
# Tested here:
#   1. the forward/inverse maps and the single-domain interval helper (bounds, hand computation);
#   2. point estimates and variances do not depend on interval_scale;
#   3. transformed bounds are valid with no clipping (arcsin in [0, 1], log positive), where the natural
#      ones sit on the bounds / go negative;
#   4. keep_replicates stores the matrices (dimensions, labels, perturbed vs unperturbed), and rebuilding the
#      interval from them with the helper reproduces the package's interval exactly, on both scales;
#   5. the argument checks.

make_scale_data <- function(seed = 21, positive = FALSE) {
  set.seed(seed)
  nd <- 24; ns <- 10
  pop <- expand.grid(sa = seq_len(ns), dom = seq_len(nd))
  pop$sub <- sprintf("%02d_%02d", pop$dom, pop$sa)
  pop$reg <- sprintf("R%d", (pop$dom - 1) %/% 8 + 1)
  pop$x1 <- rnorm(nrow(pop)); pop$x2 <- runif(nrow(pop))
  u <- rnorm(nd, 0, 0.5)[pop$dom]
  pop$npop <- round(runif(nrow(pop), 50, 500))
  ## a low rate, so that many wards sit near the lower bound
  pop$p <- plogis(-2.5 + 0.9 * pop$x1 + 0.5 * pop$x2 + u)
  pop$dom <- sprintf("%02d", pop$dom)
  sd_ <- unlist(lapply(split(unique(pop$dom), (seq_len(nd) - 1) %/% 8), function(d) sort(sample(d, 5))))
  smp <- do.call(rbind, lapply(sd_, function(d) { s <- pop[pop$dom == d, ][sample(ns, 4), ]; s }))
  smp$y <- if (positive) 0.05 + rbinom(nrow(smp), 25, smp$p) / 25 else rbinom(nrow(smp), 25, smp$p) / 25
  smp$wt <- exp(rnorm(nrow(smp), 0, 0.5)) * 100
  list(pop = pop[, c("sub", "dom", "reg", "x1", "x2", "npop")], smp = smp[, c("sub", "dom", "reg", "x1", "x2", "y", "wt")])
}

fit_scale <- function(d, transformation = "arcsin", interval_scale = "natural", keep_replicates = FALSE,
                      perturb = TRUE, B = 40, benchmark = "Mean") {
  args <- list(fixed = y ~ x1 + x2, smp_data = d$smp, smp_weights = "wt",
               pop_data = d$pop, pop_weights = "npop", domains = "dom", sub_domains = "sub",
               transformation = transformation, bootstrap = TRUE, B = B, L = 5, cpus = 1, seed = 3,
               nrounds = 20, eta = 0.3, max_depth = 2, colsample_bylevel = 1, colsample_bynode = 1, subsample = 1,
               interval_scale = interval_scale, keep_replicates = keep_replicates,
               interval_method = "plus")      # these tests are about the scale; the default method is tested in test_xgb_interval_method.R
  if (!is.null(benchmark)) {
    args <- c(args, list(benchmark = benchmark, benchmark_level = "reg",
                         benchmark_type = if (transformation == "log") "ratio" else "logit_raking"))
    if (perturb) args <- c(args, list(perturb_benchmark = TRUE,
                                      benchmark_target_se = c(R1 = 0.02, R2 = 0.02, R3 = 0.02)))
  }
  suppressWarnings(suppressMessages(do.call(povmap::xgb, args)))
}

test_that("the transformation maps are valid on the whole line and invert each other", {
  f <- povmap:::.interval_scale_fns("arcsin")
  expect_true(all(f$from(c(-5, -0.1, 0, 0.7, pi / 2, 3, 100)) >= 0 & f$from(c(-5, -0.1, 0, 0.7, pi / 2, 3, 100)) <= 1))
  p <- c(0, 0.001, 0.2, 0.5, 0.999, 1)
  expect_equal(f$from(f$to(p)), p, tolerance = 1e-12)
  expect_equal(f$to(c(-0.1, 1.1)), c(0, pi / 2))                 # out-of-range replicates are floored/capped, not NaN
  g <- povmap:::.interval_scale_fns("log")
  expect_true(all(g$from(c(-50, 0, 3)) > 0))
  expect_equal(g$from(g$to(c(0.01, 1, 40))), c(0.01, 1, 40), tolerance = 1e-12)
  expect_true(is.finite(g$to(0)))                                # a replicate of exactly 0 does not give -Inf
  s <- povmap:::.interval_scale_fns("sqrt")
  expect_equal(s$from(s$to(c(0, 0.04, 0.5))), c(0, 0.04, 0.5), tolerance = 1e-12)
  expect_true(all(s$from(c(-2, 0, 2)) >= 0))
  expect_equal(povmap:::.interval_scale_fns("no")$to(3), 3)
  expect_error(povmap:::.interval_scale_fns("logistic"), "implemented for")
})

test_that("the single-domain helper reproduces a hand computation and keeps bounds valid", {
  set.seed(1)
  fns <- povmap:::.interval_scale_fns("arcsin")
  est <- 0.02
  rep_est <- pmin(pmax(est + rnorm(400, 0, 0.03), 0), 1)         # replicates that pile up at 0
  rep_truth <- pmin(pmax(est + rnorm(400, 0, 0.03), 0), 1)
  ci <- povmap:::.interval_transformed(est, rep_est, rep_truth, fns, 0.95)
  d <- fns$to(rep_est) - fns$to(rep_truth); d <- d - mean(d)
  q <- unname(quantile(d, c(0.025, 0.975)))
  expect_equal(unname(ci), c(sin(max(asin(sqrt(est)) + q[1], 0))^2, sin(min(asin(sqrt(est)) + q[2], pi / 2))^2), tolerance = 1e-12)
  expect_true(ci[["lower"]] >= 0 && ci[["upper"]] <= 1 && ci[["lower"]] <= ci[["upper"]])
  ## unbenchmarked form: spread of the replicates around their own mean
  ci0 <- povmap:::.interval_transformed(est, rep_est, NULL, fns, 0.9)
  d0 <- fns$to(rep_est); d0 <- d0 - mean(d0); q0 <- unname(quantile(d0, c(0.05, 0.95)))
  expect_equal(unname(ci0), c(sin(max(asin(sqrt(est)) + q0[1], 0))^2, sin(min(asin(sqrt(est)) + q0[2], pi / 2))^2), tolerance = 1e-12)
  ## log: always positive, whatever the spread
  ci_log <- povmap:::.interval_transformed(0.1, exp(rnorm(300, log(0.1), 1.5)), exp(rnorm(300, log(0.1), 1.5)),
                                           povmap:::.interval_scale_fns("log"), 0.95)
  expect_true(all(ci_log > 0))
})

test_that("point estimates and variances do not depend on interval_scale", {
  d <- make_scale_data()
  nat <- fit_scale(d, "arcsin", "natural"); tr <- fit_scale(d, "arcsin", "transformed")
  expect_identical(tr$ind, nat$ind)
  expect_identical(tr$var, nat$var)
  expect_identical(tr$Var_bench_unperturbed, nat$Var_bench_unperturbed)
  expect_identical(tr$yhat, nat$yhat)
})

test_that("arcsin: transformed benchmarked intervals are valid with no clipping and are not piled on the bounds", {
  d <- make_scale_data()
  nat <- fit_scale(d, "arcsin", "natural"); tr <- fit_scale(d, "arcsin", "transformed")
  for (col in c("Lower_bench", "Upper_bench", "Lower", "Upper")) {
    expect_true(all(is.finite(tr$CI[[col]])))
    expect_true(all(tr$CI[[col]] >= 0 & tr$CI[[col]] <= 1))
  }
  expect_true(all(tr$CI$Lower_bench <= tr$CI$Upper_bench))
  expect_true(all(tr$CI$Lower <= tr$CI$Upper))
  ## natural: a low rate with a perturbed benchmark piles the lower bound on 0; transformed does not
  expect_gt(mean(nat$CI$Lower_bench == 0), 0)
  expect_lt(mean(tr$CI$Lower_bench == 0), mean(nat$CI$Lower_bench == 0))
  expect_true(all(tr$CI$Lower_bench < tr$ind$Mean_bench | tr$CI$Lower_bench == 0))
  expect_true(all(tr$CI$Upper_bench > tr$ind$Mean_bench))
})

test_that("log: transformed intervals are positive where the natural ones go negative", {
  d <- make_scale_data(seed = 5, positive = TRUE)
  d$smp$y <- d$smp$y * 0.4                                         # an index-like outcome, some wards near 0.02
  d$smp$y <- pmax(d$smp$y, 0.01)
  nat <- fit_scale(d, "log", "natural", perturb = FALSE); tr <- fit_scale(d, "log", "transformed", perturb = FALSE)
  expect_identical(tr$ind, nat$ind)
  expect_true(all(tr$CI$Lower_bench > 0)); expect_true(all(tr$CI$Lower > 0))
  expect_true(all(tr$CI$Lower_bench < tr$CI$Upper_bench))
  expect_true(all(tr$CI$Lower_bench < tr$ind$Mean_bench) && all(tr$CI$Upper_bench > tr$ind$Mean_bench))
})

test_that("keep_replicates stores the matrices and the intervals can be rebuilt from them", {
  d <- make_scale_data()
  B <- 40
  tr <- fit_scale(d, "arcsin", "transformed", keep_replicates = TRUE, B = B)
  nat <- fit_scale(d, "arcsin", "natural", keep_replicates = TRUE, B = B)
  expect_null(fit_scale(d, "arcsin", "natural", keep_replicates = FALSE, B = B)$replicates)
  r <- tr$replicates
  nd <- nrow(tr$ind)
  expect_identical(sort(r$domains), sort(as.character(tr$ind$Domain)))
  for (m in c("estimate", "benchmarked", "benchmarked_unperturbed", "truth")) {
    expect_equal(dim(r[[m]]), c(B, nd), info = m)
    expect_identical(colnames(r[[m]]), r$domains, info = m)
  }
  expect_equal(dim(r$benchmark_perturbation), c(B, 3))
  expect_gt(max(abs(r$benchmarked - r$benchmarked_unperturbed)), 0)    # the perturbation changed the replicates
  ## replicates do not depend on interval_scale
  expect_identical(nat$replicates$benchmarked, r$benchmarked)
  expect_identical(nat$replicates$truth, r$truth)
  ## rebuild both intervals from the stored replicates
  fns <- povmap:::.interval_scale_fns("arcsin")
  i <- match(as.character(tr$ind$Domain), r$domains)
  for (k in c(1, 7, nd)) {
    ci_b <- povmap:::.interval_transformed(tr$ind$Mean_bench[k], r$benchmarked[, i[k]], r$truth[, i[k]], fns, 0.95)
    expect_equal(unname(ci_b), c(tr$CI$Lower_bench[k], tr$CI$Upper_bench[k]), tolerance = 1e-12)
    ci_u <- povmap:::.interval_transformed(tr$ind$Mean[k], r$estimate[, i[k]], NULL, fns, 0.95)
    expect_equal(unname(ci_u), c(tr$CI$Lower[k], tr$CI$Upper[k]), tolerance = 1e-12)
    ## and the natural-scale interval from the same matrices
    e <- r$benchmarked[, i[k]] - r$truth[, i[k]]; e <- e - mean(e)
    nat_b <- pmin(pmax(nat$ind$Mean_bench[k] + quantile(e, c(0.025, 0.975)), 0), 1)
    expect_equal(unname(nat_b), c(nat$CI$Lower_bench[k], nat$CI$Upper_bench[k]), tolerance = 1e-12)
  }
  ## without a perturbation the two benchmarked matrices coincide; without a benchmark there are none
  np <- fit_scale(d, "arcsin", "natural", keep_replicates = TRUE, perturb = FALSE, B = B)
  expect_identical(np$replicates$benchmarked, np$replicates$benchmarked_unperturbed)
  nb <- fit_scale(d, "arcsin", "natural", keep_replicates = TRUE, benchmark = NULL, B = B)
  expect_null(nb$replicates$benchmarked); expect_null(nb$replicates$truth)
  expect_equal(dim(nb$replicates$estimate), c(B, nd))
})

test_that("interval_scale is checked", {
  d <- make_scale_data()
  expect_error(fit_scale(d, "arcsin", "bogus"), "should be one of")
  expect_error(suppressWarnings(povmap::xgb(
    fixed = y ~ x1 + x2, smp_data = d$smp, smp_weights = "wt", pop_data = d$pop, pop_weights = "npop",
    domains = "dom", sub_domains = "sub", transformation = "arcsin", B = 4, L = 2, nrounds = 5,
    interval_scale = "transformed", boot_estimates = TRUE)), "boot_estimates")
  expect_error(suppressWarnings(povmap::xgb(
    fixed = y ~ x1 + x2, smp_data = d$smp, smp_weights = "wt", pop_data = d$pop, pop_weights = "npop",
    domains = "dom", sub_domains = "sub", transformation = "logistic", B = 4, L = 2, nrounds = 5,
    interval_scale = "transformed")), "implemented for")
})
