test_that("smearing = 'naive' returns lambda 1, so the RNG stream is untouched", {
  # The caller multiplies the drawn residual by this value. Returning exactly 1
  # is what makes smearing = "naive" bit-identical to the behaviour before the
  # option existed: same sample() call, same stream, multiplied by a no-op.
  lam <- povmap:::.xgb_smearing_lambda(
    resid = rnorm(100), wt = rep(1, 100), var_y = rep(0.01, 100),
    y = rep(0.05, 100), transform_outcome = povmap:::sqrt_transform,
    smearing = "naive")
  expect_identical(lam, 1)
  # and it is a no-op even when variance_y is absent, which is the default path
  expect_identical(
    povmap:::.xgb_smearing_lambda(resid = rnorm(10), wt = rep(1, 10), var_y = NULL,
                                  y = rep(0.05, 10),
                                  transform_outcome = povmap:::sqrt_transform),
    1)
})

test_that("smearing = 'corrected' errors when variance_y is absent", {
  # Erroring rather than silently falling back: the sampling component cannot be
  # identified from the residuals alone, so a quiet default would return a
  # plausible number computed from the wrong quantity.
  expect_error(
    povmap:::.xgb_smearing_lambda(
      resid = rnorm(50), wt = rep(1, 50), var_y = NULL, y = rep(0.05, 50),
      transform_outcome = povmap:::sqrt_transform, smearing = "corrected"),
    "requires variance_y")
})

test_that("lambda equals sqrt(1 - var_sampling/var_resid) on the transformed scale", {
  # sqrt transform: dz/dy = 1/(2 sqrt(y)), so the delta-method map is exact enough
  # to check the arithmetic against a hand computation.
  set.seed(4)
  n <- 20000; y <- rep(0.04, n); wt <- rep(1, n)
  vy <- rep(0.0004, n)                       # rate-scale sampling variance
  dzdy <- 1 / (2 * sqrt(0.04))               # = 2.5
  var_samp_t <- 0.0004 * dzdy^2              # = 0.0025
  resid <- rnorm(n, 0, 0.10)                 # transformed-scale residual sd 0.10
  lam <- povmap:::.xgb_smearing_lambda(
    resid = resid, wt = wt, var_y = vy, y = y,
    transform_outcome = povmap:::sqrt_transform, smearing = "corrected")
  expect_equal(lam, sqrt(max(0, 1 - var_samp_t / mean(resid^2))), tolerance = 1e-8)
  expect_gt(lam, 0); expect_lt(lam, 1)
})

test_that("lambda goes to 0 when sampling variance dominates, falling back to g(z)", {
  # The degrade-safely property. If the residual pool is all sampling noise there
  # is no true between-cell variation to smear over, and the estimator should
  # collapse to the back-transform of the mean rather than return something
  # arbitrary or negative under the square root.
  set.seed(5)
  n <- 5000
  lam <- povmap:::.xgb_smearing_lambda(
    resid = rnorm(n, 0, 0.01), wt = rep(1, n),
    var_y = rep(1, n),                        # sampling variance far exceeding it
    y = rep(0.04, n), transform_outcome = povmap:::sqrt_transform,
    smearing = "corrected")
  expect_identical(lam, 0)
  expect_false(is.na(lam))
})

test_that("the correction removes the sqrt smearing bias it is meant to remove", {
  # End-to-end on the arithmetic the correction targets: with residual variance
  # made up of a true part and a sampling part, smearing the FULL pool overshoots
  # z^2 + var_true, and scaling by lambda recovers it.
  set.seed(6)
  n <- 2e6
  z <- sqrt(0.006)
  sd_true <- 0.04; sd_samp <- 0.03
  sd_tot  <- sqrt(sd_true^2 + sd_samp^2)
  target  <- z^2 + sd_true^2                  # what smearing SHOULD return
  e <- rnorm(n, 0, sd_tot)
  naive <- mean(povmap:::sqrt_transform_back(z + e))
  lam   <- sqrt(max(0, 1 - sd_samp^2 / sd_tot^2))
  corrected <- mean(povmap:::sqrt_transform_back(z + lam * e))
  expect_gt(naive, target)                                  # naive overshoots
  expect_equal(corrected, target, tolerance = 2e-4)         # corrected lands on it
  expect_lt(abs(corrected - target), abs(naive - target))
})
