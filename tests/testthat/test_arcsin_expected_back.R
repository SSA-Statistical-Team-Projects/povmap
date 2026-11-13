test_that("arcsin_expected_back matches numerical integration of sin(x)^2 under normal error", {
  mus <- seq(0.02, pi / 2 - 0.02, length.out = 15)
  sds <- c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5)
  for (s in sds) {
    num <- vapply(mus, function(m) {
      integrate(function(e) sin(m + e)^2 * dnorm(e, 0, s), -Inf, Inf, rel.tol = 1e-10)$value
    }, numeric(1))
    expect_equal(povmap:::arcsin_expected_back(mus, s^2), num, tolerance = 1e-8)
  }
})

test_that("arcsin_expected_back reduces to sin(mu)^2 at zero variance and to 1/2 at large variance", {
  mu <- c(0.1, 0.7, 1.4)
  expect_equal(povmap:::arcsin_expected_back(mu, 0), sin(mu)^2, tolerance = 1e-14)
  expect_equal(povmap:::arcsin_expected_back(mu, 50), rep(0.5, 3), tolerance = 1e-12)
})

test_that("expected_untransformed_mean uses the exact arcsin form", {
  mu <- c(0.3, 0.9, 1.3); v <- 0.04
  expect_equal(povmap:::expected_untransformed_mean(mu = mu, var = v, transformation = "arcsin", lambda = NULL),
               0.5 - 0.5 * cos(2 * mu) * exp(-2 * v), tolerance = 1e-14)
})

test_that("arcsin_expected_back clamps mu to [0, pi/2] and is monotone in mu", {
  mu <- seq(-1, 3, length.out = 401)
  for (v in c(0, 0.01, 0.1, 0.25)) {
    e <- povmap:::arcsin_expected_back(mu, v)
    expect_true(all(diff(e) >= -1e-14))
    expect_equal(e[mu <= 0], rep(povmap:::arcsin_expected_back(0, v), sum(mu <= 0)))
    expect_equal(e[mu >= pi / 2], rep(povmap:::arcsin_expected_back(pi / 2, v), sum(mu >= pi / 2)))
  }
  expect_equal(povmap:::arcsin_expected_back(c(-0.5, 2), 0), c(0, 1))
})
