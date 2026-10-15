sim_megb_lmm <- function(weighted) {
  set.seed(11)
  D <- 12
  n_i <- c(1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 15, 20)
  d <- data.frame(dom = factor(rep(sprintf("d%02d", seq_len(D)), n_i)))
  d$r <- rnorm(D, 0, 0.6)[as.integer(d$dom)] + rnorm(nrow(d), 0, 0.5)
  d$w <- if (weighted) runif(nrow(d), 0.3, 3) else 1
  list(d = d, fit = suppressMessages(lme4::lmer(r ~ 1 + (1 | dom), data = d, REML = FALSE,
                                                weights = if (weighted) d$w else NULL)))
}

test_that(".megb_re_condvar equals sigma_b^2 (1 - gamma_i) for sampled and sigma_b^2 for unsampled domains", {
  for (weighted in c(FALSE, TRUE)) {
    s <- sim_megb_lmm(weighted)
    sb2 <- as.data.frame(lme4::VarCorr(s$fit))$vcov[1]; se2 <- stats::sigma(s$fit)^2
    W <- tapply(s$d$w, s$d$dom, sum)                       # n_i when unweighted
    gam <- sb2 * W / (sb2 * W + se2)
    doms <- c(levels(s$d$dom), "new1", "new2")
    v <- povmap:::.megb_re_condvar(s$fit, doms)
    expect_equal(v[seq_along(W)], as.numeric(sb2 * (1 - gam)), tolerance = 1e-6)
    expect_equal(v[length(W) + 1:2], rep(sb2, 2))
    expect_true(all(v[seq_along(W)] < sb2))
  }
})

# E over b ~ N(0, v_re) and e ~ N(0, sigma2_e) of g^{-1}(mu + b + e), by nested numerical integration
num_back <- function(mu, v_re, se2, g_inv) {
  inner <- function(b) vapply(b, function(bb)
    integrate(function(e) g_inv(mu + bb + e) * dnorm(e, 0, sqrt(se2)), -10 * sqrt(se2), 10 * sqrt(se2), rel.tol = 1e-10)$value,
    numeric(1))
  integrate(function(b) inner(b) * dnorm(b, 0, sqrt(v_re)), -10 * sqrt(v_re), 10 * sqrt(v_re), rel.tol = 1e-9)$value
}

test_that("corrected_bt with the random-effect variance matches numerical integration (log and arcsin)", {
  for (weighted in c(FALSE, TRUE)) {
    s <- sim_megb_lmm(weighted)
    sb2 <- as.data.frame(lme4::VarCorr(s$fit))$vcov[1]; se2 <- stats::sigma(s$fit)^2
    v <- povmap:::.megb_re_condvar(s$fit, c(levels(s$d$dom), "unsampled"))
    pick <- c(1, 6, 12, 13)                                  # small, medium, large sampled domain; unsampled
    for (tr in c("log", "arcsin")) {
      cbt <- povmap:::.megb_make_corrected_bt(se2, tr)
      g_inv <- if (tr == "log") exp else function(x) sin(x)^2
      mu <- if (tr == "log") c(-0.4, 0.3, 1.1, 0.7) else c(0.3, 0.8, 1.2, 1.5)
      for (k in seq_along(pick)) {
        expect_equal(cbt(mu[k], v[pick[k]]), num_back(mu[k], v[pick[k]], se2, g_inv), tolerance = 1e-7)
      }
      # the sampled-domain correction is smaller than the unsampled one, and both exceed the unit-error-only form
      expect_lt(cbt(mu[1], v[1]), cbt(mu[1], v[13]))
    }
  }
})

test_that("re_var = 0 reproduces the unit-error-only back-transformation", {
  expect_equal(povmap:::.megb_make_corrected_bt(0.2, "log")(0.5), exp(0.5 + 0.1))
  expect_equal(povmap:::.megb_make_corrected_bt(0.2, "arcsin")(0.5), povmap:::arcsin_expected_back(0.5, 0.2))
  expect_equal(povmap:::.megb_make_corrected_bt(0, "log")(0.5, 0.3), exp(0.5 + 0.15))
})

test_that("megb runs with the random-effect variance correction, point estimate and bootstrap, log and arcsin", {
  set.seed(7); D <- 14; n_pop <- 15; n_smp <- 6
  pop <- data.frame(dom = factor(rep(sprintf("d%02d", 1:D), each = n_pop)), x1 = rnorm(D * n_pop), x2 = runif(D * n_pop), w = rpois(D * n_pop, 20) + 1)
  u <- rnorm(D, 0, 0.4)[as.integer(pop$dom)]
  pop$y <- exp(0.4 * pop$x1 + u + rnorm(nrow(pop), 0, 0.4))
  pop$p <- plogis(0.5 * pop$x1 + u + rnorm(nrow(pop), 0, 0.3))
  smp <- do.call(rbind, lapply(split(pop, pop$dom)[1:11], function(d) d[sample(nrow(d), n_smp), ]))
  gp <- list(eta = 0.3, max_depth = 2, subsample = 1, nrounds = 30)
  for (cfg in list(list(f = y ~ x1, tr = "log"), list(f = p ~ x1, tr = "arcsin"))) {
    m <- suppressMessages(suppressWarnings(megb(
      fixed = update(cfg$f, . ~ x1 + x2), smp_data = smp, pop_data = pop[, c("dom", "x1", "x2", "w")], pop_weights = "w",
      domains = "dom", transformation = cfg$tr, mse = TRUE, B = 3, bootstrap_refit = "lmm_only",
      na.rm = FALSE, seed = 5, gradient_params = gp)))
    expect_false(anyNA(m$ind$Mean))
    expect_false(anyNA(m$MSE$Mean %||% m$var$Mean))
  }
})

test_that("the correction varies with n_i and matches numerical integration for every domain (n_i from 2 to 50)", {
  set.seed(21)
  D <- 15
  n_i <- round(seq(2, 50, length.out = D))
  d <- data.frame(dom = factor(rep(sprintf("d%02d", seq_len(D)), n_i)))
  d$r <- rnorm(D, 0, 0.6)[as.integer(d$dom)] + rnorm(nrow(d), 0, 0.5)
  for (weighted in c(FALSE, TRUE)) {
    d$w <- if (weighted) runif(nrow(d), 0.3, 3) else 1
    fit <- suppressMessages(lme4::lmer(r ~ 1 + (1 | dom), data = d, REML = FALSE,
                                       weights = if (weighted) d$w else NULL))
    se2 <- stats::sigma(fit)^2
    v <- povmap:::.megb_re_condvar(fit, levels(d$dom))
    expect_true(all(diff(v) < 0 | weighted))                  # unweighted: variance falls as n_i grows
    expect_gt(diff(range(v)), 0.05 * max(v))                  # and it does vary across domains
    for (tr in c("log", "arcsin")) {
      cbt <- povmap:::.megb_make_corrected_bt(se2, tr)
      g_inv <- if (tr == "log") exp else function(x) sin(x)^2
      mu <- if (tr == "log") 0.4 else 0.9
      num <- vapply(seq_len(D), function(k) num_back(mu, v[k], se2, g_inv), numeric(1))
      expect_equal(cbt(rep(mu, D), v), num, tolerance = 1e-7)
    }
  }
})

test_that(".megb_re_var computes the random-effect variance only for log, log.shift, poisson and arcsin", {
  s <- sim_megb_lmm(FALSE)
  doms <- levels(s$d$dom)
  for (tr in c("log", "log.shift", "poisson", "arcsin")) {
    cbt <- povmap:::.megb_make_corrected_bt(0.2, tr, if (tr == "log.shift") 1 else NULL)
    expect_equal(povmap:::.megb_re_var(cbt, s$fit, doms), povmap:::.megb_re_condvar(s$fit, doms))
  }
  for (tr in c("no", "logistic")) {
    cbt <- povmap:::.megb_make_corrected_bt(0.2, tr)
    expect_identical(povmap:::.megb_re_var(cbt, s$fit, doms), rep(0, length(doms)))
  }
})

# Untransformed and logistic megb are unchanged from the v0.6 convention release (b64e60a, povmap 1.0.2,
# xgboost 3.1.2.1). The fixture holds that build's ind, var (the MSE-type uncertainty) and CI on this data.
test_that("megb(transformation = 'no' / 'logistic') matches the b64e60a fixture", {
  fx <- testthat::test_path("fixtures", "megb_identity_b64e60a.rds")
  skip_if_not(file.exists(fx), "fixture missing")
  ref <- readRDS(fx)
  set.seed(7); D <- 14; n_pop <- 15; n_smp <- 6
  pop <- data.frame(dom = factor(rep(sprintf("d%02d", 1:D), each = n_pop)), x1 = rnorm(D * n_pop), x2 = runif(D * n_pop), w = rpois(D * n_pop, 20) + 1)
  pop$y <- 0.4 * pop$x1 + rnorm(D, 0, 0.4)[as.integer(pop$dom)] + rnorm(nrow(pop), 0, 0.4)
  pop$p <- plogis(0.5 * pop$x1 + rnorm(D, 0, 0.4)[as.integer(pop$dom)] + rnorm(nrow(pop), 0, 0.3))
  smp <- do.call(rbind, lapply(split(pop, pop$dom)[1:11], function(d) d[sample(nrow(d), n_smp), ]))
  for (cfg in list(list(f = y ~ x1 + x2, tr = "no"), list(f = p ~ x1 + x2, tr = "logistic"))) {
    m <- suppressMessages(suppressWarnings(megb(fixed = cfg$f, smp_data = smp, pop_data = pop[, c("dom", "x1", "x2", "w")],
      pop_weights = "w", domains = "dom", transformation = cfg$tr, mse = TRUE, B = 3, bootstrap_refit = "lmm_only",
      na.rm = FALSE, seed = 5, gradient_params = list(eta = 0.3, max_depth = 2, subsample = 1, nrounds = 30))))
    expect_identical(m$ind, ref[[cfg$tr]]$ind)
    expect_equal(m$var, ref[[cfg$tr]]$var, tolerance = 1e-10)
    expect_equal(m$CI, ref[[cfg$tr]]$CI, tolerance = 1e-10)
  }
})
