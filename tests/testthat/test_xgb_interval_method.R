# Tests for xgb(interval_method = "basic").
#
# The benchmarked interval is built from the centred bootstrap error e = (benchmarked replicate) - (simulated truth).
# "plus" (the default before this change) is estimate + q(alpha/2) .. estimate + q(1 - alpha/2). With e defined as estimate minus truth the
# basic bootstrap interval for the truth is estimate - q(1 - alpha/2) .. estimate - q(alpha/2): the reflected one. The two
# agree when e is symmetric about 0; with skewed errors (near a bound) "plus" puts the long tail on the wrong side.
# "basic" is the default; interval_method = "plus" reproduces the earlier results.
#
# Tested here:
#   1. the default is "basic" and naming it changes nothing; "plus" changes only the benchmarked interval (point
#      estimates, variances and the unbenchmarked interval are identical);
#   2. basic is the reflection of plus around the estimate for the quantiles used: lower + upper of "basic" and "plus" are
#      related as hand-computed from the stored replicates, on the natural and the transformed scale;
#   3. the helper with method = "basic", and the argument check.

make_method_data <- function(seed = 41) {
  set.seed(seed)
  nd <- 24; ns <- 8
  pop <- expand.grid(sa = seq_len(ns), dom = seq_len(nd))
  pop$sub <- sprintf("%02d_%02d", pop$dom, pop$sa)
  pop$reg <- sprintf("R%d", (pop$dom - 1) %/% 8 + 1)
  pop$x1 <- rnorm(nrow(pop)); pop$x2 <- runif(nrow(pop))
  u <- rnorm(nd, 0, 0.5)[pop$dom]
  pop$npop <- round(runif(nrow(pop), 50, 500))
  pop$p <- plogis(-2.5 + 0.9 * pop$x1 + 0.5 * pop$x2 + u)          # low rates: skewed errors near 0
  pop$dom <- sprintf("%02d", pop$dom)
  sd_ <- unlist(lapply(split(unique(pop$dom), (seq_len(nd) - 1) %/% 8), function(d) sort(sample(d, 5))))
  smp <- do.call(rbind, lapply(sd_, function(d) pop[pop$dom == d, ][sample(ns, 4), ]))
  smp$y <- rbinom(nrow(smp), 25, smp$p) / 25
  smp$wt <- exp(rnorm(nrow(smp), 0, 0.5)) * 100
  list(pop = pop[, c("sub", "dom", "reg", "x1", "x2", "npop")], smp = smp[, c("sub", "dom", "reg", "x1", "x2", "y", "wt")])
}
fit_method <- function(d, interval_method = NULL, interval_scale = "natural", B = 60, ...) {   # NULL: the package default
  extra <- if (is.null(interval_method)) list() else list(interval_method = interval_method)
  suppressWarnings(suppressMessages(do.call(povmap::xgb, c(list(
    fixed = y ~ x1 + x2, smp_data = d$smp, smp_weights = "wt", pop_data = d$pop, pop_weights = "npop", domains = "dom", sub_domains = "sub",
    transformation = "arcsin", bootstrap = TRUE, B = B, L = 5, cpus = 1, seed = 3, nrounds = 20, eta = 0.3, max_depth = 2,
    colsample_bylevel = 1, colsample_bynode = 1, subsample = 1, benchmark = "Mean", benchmark_level = "reg", benchmark_type = "logit_raking",
    perturb_benchmark = TRUE, benchmark_target_se = c(R1 = 0.01, R2 = 0.01, R3 = 0.01), keep_replicates = TRUE,
    interval_scale = interval_scale), extra, list(...)))))
}

test_that("the default is basic, and plus (the earlier default) changes only the benchmarked interval", {
  d <- make_method_data()
  a <- fit_method(d); a2 <- fit_method(d, "basic"); p <- fit_method(d, "plus")
  expect_identical(a$CI, a2$CI); expect_identical(a$ind, a2$ind)                           # naming the default changes nothing
  expect_identical(p$ind, a$ind); expect_identical(p$var, a$var)
  expect_identical(p$CI$Lower, a$CI$Lower); expect_identical(p$CI$Upper, a$CI$Upper)      # unbenchmarked interval
  expect_false(isTRUE(all.equal(p$CI$Lower_bench, a$CI$Lower_bench)))
  expect_true(all(a$CI$Lower_bench >= 0 & a$CI$Upper_bench <= 1))
  expect_true(all(p$CI$Lower_bench >= 0 & p$CI$Upper_bench <= 1))
})

test_that("basic is the reflected interval, rebuilt by hand from the stored replicates", {
  d <- make_method_data()
  for (sc in c("natural", "transformed")) {
    fb <- fit_method(d, "basic", sc); fp <- fit_method(d, "plus", sc)
    r <- fb$replicates; i <- match(as.character(fb$ind$Domain), r$domains)
    to <- if (sc == "natural") identity else function(y) asin(sqrt(pmin(pmax(y, 0), 1)))
    from <- if (sc == "natural") function(t) pmin(pmax(t, 0), 1) else function(t) sin(pmax(0, pmin(t, pi / 2)))^2
    for (k in c(1, 9, nrow(fb$ind))) {
      e <- to(r$benchmarked[, i[k]]) - to(r$truth[, i[k]]); e <- e - mean(e)
      q <- unname(quantile(e, c(0.025, 0.975))); th <- to(fb$ind$Mean_bench[k])
      expect_equal(c(fb$CI$Lower_bench[k], fb$CI$Upper_bench[k]), from(c(th - q[2], th - q[1])), tolerance = 1e-12, info = sc)
      expect_equal(c(fp$CI$Lower_bench[k], fp$CI$Upper_bench[k]), from(c(th + q[1], th + q[2])), tolerance = 1e-12, info = sc)
    }
    ## with skewed errors the two differ in where the long tail falls
    expect_gt(max(abs((fb$CI$Upper_bench - fb$ind$Mean_bench) - (fp$CI$Upper_bench - fp$ind$Mean_bench))), 1e-4)
  }
})

test_that("the helper and the argument check", {
  fns <- povmap:::.interval_scale_fns("no")
  set.seed(2)
  est <- 1; rep_est <- 1 + rexp(500, 3); rep_truth <- 1 + rexp(500, 3)            # skewed errors
  e <- rep_est - rep_truth; e <- e - mean(e); q <- unname(quantile(e, c(0.025, 0.975)))
  expect_equal(unname(povmap:::.interval_transformed(est, rep_est, rep_truth, fns, 0.95, "plus")), c(est + q[1], est + q[2]), tolerance = 1e-12)
  expect_equal(unname(povmap:::.interval_transformed(est, rep_est, rep_truth, fns, 0.95, "basic")), c(est - q[2], est - q[1]), tolerance = 1e-12)
  ## symmetric errors: the two coincide
  s <- c(-2, -1, 0, 1, 2); expect_equal(unname(povmap:::.interval_transformed(0, s, rep(0, 5), fns, 0.8, "plus")),
                                         unname(povmap:::.interval_transformed(0, s, rep(0, 5), fns, 0.8, "basic")))
  d <- make_method_data()
  expect_error(fit_method(d, "bogus"), "should be one of")
})
