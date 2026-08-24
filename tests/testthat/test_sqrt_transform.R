test_that("sqrt_transform is defined at zero and round-trips exactly", {
  y <- c(0, 1e-12, 1e-8, 0.005, 0.05, 0.5, 0.99, 1)
  z <- povmap:::sqrt_transform(y)$y
  expect_true(all(is.finite(z)))
  expect_identical(z[1], 0)                      # defined AT zero, no epsilon needed
  back <- povmap:::sqrt_transform_back(z)
  expect_equal(back, y, tolerance = 1e-14)
  # the case arcsin can only handle by clamping, sqrt handles natively
  expect_identical(povmap:::sqrt_transform_back(povmap:::sqrt_transform(0)$y), 0)
})

test_that("sqrt_transform_back is WEAKLY MONOTONE over the whole real line", {
  # This is the property the clamp exists for. z^2 decreases for z < 0, so an
  # unclamped inverse would fold back upward and a residual draw pushing z
  # negative would RAISE the predicted rate. Same role as the [0, pi/2] clamp
  # in arcsin_transform_back.
  z <- seq(-3, 3, by = 0.001)
  b <- povmap:::sqrt_transform_back(z)
  expect_true(all(diff(b) >= -1e-15))
  expect_false(is.unsorted(b))
})

test_that("sqrt_transform_back is non-negative and clamps below zero", {
  expect_identical(povmap:::sqrt_transform_back(-1), 0)
  expect_identical(povmap:::sqrt_transform_back(-5), 0)
  expect_true(all(povmap:::sqrt_transform_back(seq(-10, 10, by = 0.01)) >= 0))
  # an unclamped inverse would give 1 and 25 here; the clamp is doing real work
  expect_false(isTRUE(all.equal(povmap:::sqrt_transform_back(-1), 1)))
})

test_that("sqrt routes through both dispatchers", {
  y <- c(0, 0.005, 0.25, 1)
  fwd <- povmap:::transformation(y = y, transformation = "sqrt", lambda = NULL,
                                 shift = NULL, framework = NULL, fixed = NULL)
  expect_equal(fwd$y, sqrt(y), tolerance = 1e-14)
  bk <- povmap:::back_transformation(y = fwd$y, transformation = "sqrt", lambda = NULL,
                                     shift = NULL, framework = NULL, fixed = NULL)
  expect_equal(bk, y, tolerance = 1e-14)
})

test_that("near zero, sqrt and arcsin are the same transformation to first order", {
  # Recorded as a check on a PREDICTION made before running the arm: because
  # asin(u) = u + u^3/6 + ..., at rates well under 1% the sqrt and arcsin scales
  # agree to several decimals, so the two arms should produce near-identical
  # estimates. If this test ever fails the prediction was wrong, not the code.
  y <- c(0.001, 0.006, 0.015)
  s <- povmap:::sqrt_transform(y)$y
  a <- povmap:::arcsin_transform(y)$y
  expect_true(all(abs(a - s) / pmax(s, 1e-12) < 0.01))   # within 1% on the transformed scale
  # and the divergence grows with y, as the series says it must
  expect_gt(abs(povmap:::arcsin_transform(0.5)$y - povmap:::sqrt_transform(0.5)$y),
            abs(povmap:::arcsin_transform(0.006)$y - povmap:::sqrt_transform(0.006)$y))
})

test_that("the smearing bias of the sqrt back-transform is exactly the residual variance", {
  # E[(z+e)^2] = z^2 + Var(e), EXACTLY -- not a second-order approximation. This
  # is what makes the correction a single subtraction, and it is the reason
  # smearing = "corrected" is worth having for this transform in particular.
  set.seed(11)
  z <- sqrt(0.006); sigma <- 0.05
  e <- rnorm(4e6, 0, sigma)
  naive <- mean(povmap:::sqrt_transform_back(z + e))
  expect_equal(naive, z^2 + sigma^2, tolerance = 5e-5)
  expect_gt(naive, z^2)          # UPWARD, because the inverse is convex
})
