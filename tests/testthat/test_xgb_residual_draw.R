# Regression tests for the sub-area residual draw in point_estim_xgb().
#
# Sub-area residuals are built by subtracting the WEIGHT-weighted domain mean,
# so sum(resid * w) == 0 holds by construction. Re-drawing those residuals is
# therefore unbiased only when the draw probability is proportional to w.
#
# R/xgb.R:1169 previously drew with prob = 1/w, the inverse of the survey
# weight. That injects a non-zero mean into every domain: in the Colombia
# migration SAE study it added +0.0127 to every domain estimate against an
# observed excess of +0.0129. The adjacent area-level draw at :1175 already used
# plain weights, and the bootstrap path at :610 already used pop_subarea_d, so
# the inverse draw disagreed both with its neighbour and with the bootstrap.

test_that("weight-centred residuals have zero weighted mean by construction", {
  set.seed(11)
  n <- 500
  w <- exp(rnorm(n, 0, 0.8))                 # unequal survey weights
  y <- 0.3 + 0.5 * runif(n) + rnorm(n, 0, 0.1)
  resid <- y - weighted.mean(y, w)           # the construction used in xgb
  expect_equal(weighted.mean(resid, w), 0, tolerance = 1e-12)
  # and the UNWEIGHTED mean is NOT zero, which is why the draw probability matters
  expect_gt(abs(mean(resid)), 1e-4)
})

test_that("drawing residuals with prob = w is unbiased, prob = 1/w is not", {
  set.seed(11)
  n <- 500
  w <- exp(rnorm(n, 0, 0.8))
  y <- 0.3 + 0.5 * runif(n) + rnorm(n, 0, 0.1)
  resid <- y - weighted.mean(y, w)

  # Expected value of a draw is sum(resid * p) with p the normalised prob.
  exp_w   <- sum(resid * (w / sum(w)))       # prob = w      -> zero
  exp_inv <- sum(resid * ((1 / w) / sum(1 / w)))  # prob = 1/w -> biased
  expect_equal(exp_w, 0, tolerance = 1e-12)
  expect_gt(abs(exp_inv), 1e-3)

  # and the same thing by simulation, so the test fails on the drawn mean itself
  set.seed(12)
  N <- 200000L
  draw_w   <- sample(resid, N, replace = TRUE, prob = w)
  draw_inv <- sample(resid, N, replace = TRUE, prob = 1 / w)
  se <- sd(resid) / sqrt(N)
  expect_lt(abs(mean(draw_w)), 4 * se)       # consistent with zero
  expect_gt(abs(mean(draw_inv)), 4 * se)     # demonstrably not zero
})

test_that("point_estim_xgb draws sub-area residuals with prob = w, not 1/w", {
  skip_if_not_installed("povmap")
  b <- deparse(povmap:::point_estim_xgb)
  hits <- grep("prob = ", b, value = TRUE)
  expect_true(length(hits) >= 1L)
  # the sub-area draw must use the weight itself
  expect_true(any(grepl("prob = sample[, fwk$smp_weights]", hits, fixed = TRUE)))
  # and must never use its inverse
  expect_false(any(grepl("prob = 1/sample", hits, fixed = TRUE)))
})
