# megb's bootstrap (5 Oct 2026):
#  * the interval subtracts the bootstrap error quantiles from the point estimate
#    (basic bootstrap); it used to add them, which is wrong when the errors are skewed;
#  * leaves_only replicates refresh the leaves once on y* by default (a documented
#    approximation); with options(povmap.leaves_only.em = TRUE) each replicate runs
#    the point fit's EM with the trees fixed: the fold models' leaves are refreshed
#    on y* net of the estimated random effect, the random effect is re-estimated
#    from the out-of-fold residuals, to the point fit's iteration cap and tolerance.
.ie_data <- function(D = 30, n_pop = 20, n_smp = 8, seed = 7) {
  set.seed(seed)
  pop <- data.frame(dom = factor(rep(sprintf("d%02d", 1:D), each = n_pop)),
                    x1 = rnorm(D * n_pop), x2 = runif(D * n_pop), w = rpois(D * n_pop, 20) + 1)
  pop$xm <- rnorm(D)[as.integer(pop$dom)]
  pop$y <- 0.5 * pop$x1 + 0.4 * pop$xm + rnorm(D, 0, 0.3)[as.integer(pop$dom)] + rnorm(nrow(pop))
  smp <- do.call(rbind, lapply(split(pop, pop$dom)[1:25], function(d) d[sample(nrow(d), n_smp), ]))
  smp$wt <- runif(nrow(smp), 0.3, 3)
  list(pop = pop, smp = smp)
}
.ie_fit <- function(M, mse = TRUE, B = 3, ...)
  suppressMessages(suppressWarnings(povmap::megb(
  fixed = y ~ x1 + x2 + xm, smp_data = M$smp, smp_weights = "wt",
  pop_data = M$pop[, c("dom", "x1", "x2", "xm", "w")], pop_weights = "w", domains = "dom",
  transformation = "no", mse = mse, B = B, na.rm = FALSE, seed = 5, cv_nfold = 5,
  gradient_params = list(eta = 0.3, max_depth = 3, subsample = 1, nrounds = 100),
  predict_sampled = "crossfit", predict_unsampled = "foldmean", bootstrap_refit = "leaves_only", ...)))

test_that("the interval subtracts the error quantiles (skewed errors)", {
  skip_if_not_installed("xgboost")
  M <- .ie_data()
  ## replace each replicate's estimate by its truth plus a right-skewed error, so the error law is known
  set.seed(11)
  orig <- utils::getFromNamespace("mse_megb", "povmap")
  skew <- function(...) {
    r <- orig(...); D <- nrow(r$tau_star_orig_unb)
    E <- matrix(stats::rexp(D * ncol(r$tau_star_orig_unb), rate = 20), D)   # mean 0.05, long right tail
    r$tau_b_orig <- r$tau_star_orig_unb + E; r$E <- E; assign("E", E, envir = rec); r
  }
  rec <- new.env()
  local_mocked_bindings(mse_megb = skew, .package = "povmap")
  f <- .ie_fit(M, B = 40)
  E <- rec$E; Ec <- E - rowMeans(E)
  q_lo <- apply(Ec, 1, stats::quantile, probs = 0.025, names = FALSE)
  q_hi <- apply(Ec, 1, stats::quantile, probs = 0.975, names = FALSE)
  ci <- merge(f$CI, f$ind, by = "Domain")
  ## domains come back in the bootstrap's (sorted) order
  o  <- match(ci$Domain, sort(as.character(unique(M$pop$dom))))
  expect_equal(ci$Lower, ci$Mean - q_hi[o], tolerance = 1e-12)
  expect_equal(ci$Upper, ci$Mean - q_lo[o], tolerance = 1e-12)
  ## the old orientation differs when the errors are skewed: the long tail is now below the estimate
  expect_true(all(ci$Mean - ci$Lower > ci$Upper - ci$Mean))
  ## and the variance is the centred errors' variance
  expect_equal(merge(f$var, f$ind, by = "Domain")$Mean.x, apply(Ec, 1, stats::var)[o], tolerance = 1e-12)
})

test_that("the error-interval helper has the basic bootstrap's coverage under skew", {
  set.seed(3); alpha <- 0.05; R <- 4000
  ## theta_hat = theta + e with skewed e; the bootstrap sees errors from the same law
  e    <- stats::rexp(R, 1) - 1
  boot <- matrix(stats::rexp(R * 199, 1) - 1, R)
  I    <- povmap:::.megb_error_interval(boot, alpha)
  th   <- 0; hat <- th + e
  cover_new <- mean(hat - I$q_hi <= th & th <= hat - I$q_lo)
  cover_old <- mean(hat + I$q_lo <= th & th <= hat + I$q_hi)
  expect_gt(cover_new, 0.93); expect_lt(cover_new, 0.97)   # nominal 0.95
  expect_lt(cover_old, 0.90)                                   # about 0.86 for this law
})

test_that("options(povmap.leaves_only.em = TRUE): replicates run the point fit's EM with the trees fixed", {
  skip_if_not_installed("xgboost")
  withr::local_options(povmap.leaves_only.em = TRUE)
  M <- .ie_data(); rec <- new.env(); rec$steps <- list(); rec$em <- list()
  o_step <- utils::getFromNamespace(".megb_fold_boot_step", "povmap")
  o_em   <- utils::getFromNamespace(".megb_leaf_em", "povmap")
  local_mocked_bindings(
    .megb_fold_boot_step = function(CFB, X_smp, X_pop, y, ...) { rec$steps[[length(rec$steps) + 1L]] <- y; o_step(CFB, X_smp, X_pop, y, ...) },
    .megb_leaf_em = function(CFB, X_smp, X_pop, y, w_ref, fit_data, lmm_formula, w_lmm, ...) {
      a <- list(...); n0 <- length(rec$steps); r <- o_em(CFB, X_smp, X_pop, y, w_ref, fit_data, lmm_formula, w_lmm, ...)
      rec$em[[length(rec$em) + 1L]] <- list(y = y, dom = as.character(fit_data[[all.vars(lmm_formula)[2]]]),
                                            max_it = a$max_iterations, tol = a$error_tolerance,
                                            targets = rec$steps[(n0 + 1L):length(rec$steps)], it = r$iterations)
      r },
    .package = "povmap")
  f <- .ie_fit(M, B = 3)
  expect_length(rec$em, 3L)
  it <- f$boot_diag$boot_em_iterations
  expect_equal(it, vapply(rec$em, `[[`, 0L, "it"))
  expect_true(all(it >= 1L & it <= 25L))
  for (e in rec$em) {
    expect_equal(e$max_it, povmap:::.megb_em_max_iterations); expect_equal(e$tol, povmap:::.megb_em_error_tolerance)
    expect_length(e$targets, e$it)                                 # one fold refresh per iteration
    for (t in e$targets) {
      d <- e$y - t                                                 # y* minus the target: the random effect,
      expect_gt(max(abs(d)), 1e-6)                                 # not zero (the target is not y* itself)
      expect_lt(max(tapply(d, e$dom, function(z) diff(range(z)))), 1e-8)   # and constant within a domain
    }
  }
})

test_that("the default is the single pass on y*; the EM leaves the point estimates unchanged", {
  skip_if_not_installed("xgboost")
  M <- .ie_data(); rec <- new.env(); rec$steps <- list()
  o_step <- utils::getFromNamespace(".megb_fold_boot_step", "povmap")
  local_mocked_bindings(.megb_fold_boot_step = function(CFB, X_smp, X_pop, y, ...) { rec$steps[[length(rec$steps) + 1L]] <- y; o_step(CFB, X_smp, X_pop, y, ...) },
                        .package = "povmap")
  f <- .ie_fit(M, B = 2)                                           # default: no option set
  expect_equal(f$boot_diag$boot_em_iterations, c(0L, 0L))
  expect_length(rec$steps, 2L)
  ## and the point estimate is the same with or without the EM in the bootstrap
  withr::local_options(povmap.leaves_only.em = TRUE)
  expect_identical(.ie_fit(M, B = 2)$ind, f$ind)
})
