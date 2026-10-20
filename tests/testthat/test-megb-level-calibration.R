# boot_synthetic and megb_level_calibration() (MEGB_BOOT_FIX_PREREG_v8_final.md, arm D, option A).
.lc_data <- function(D = 24, n_pop = 20, n_smp = 6, seed = 9) {
  set.seed(seed)
  pop <- data.frame(dom = factor(sprintf("d%02d", rep(seq_len(D), each = n_pop))), x1 = rnorm(D * n_pop), x2 = runif(D * n_pop), w = rpois(D * n_pop, 20) + 1)
  pop$y <- exp(0.3 * pop$x1 + rnorm(D, 0, 0.3)[as.integer(pop$dom)] + rnorm(nrow(pop), 0, 0.5))
  keep <- unique(as.character(pop$dom))[1:18]
  smp <- do.call(rbind, lapply(keep, function(d) { z <- pop[pop$dom == d, ]; z[sample(nrow(z), n_smp), ] }))
  smp$wt <- runif(nrow(smp), 0.5, 2); rownames(smp) <- NULL
  list(pop = pop[, c("dom", "x1", "x2", "w")], smp = smp)
}
.lc_fit <- function(M, ..., tr = "no", pw = NULL, synth = TRUE)
  suppressMessages(suppressWarnings(povmap::megb(fixed = y ~ x1 + x2, smp_data = M$smp, smp_weights = "wt", pop_data = M$pop, pop_weights = pw, domains = "dom",
    transformation = tr, mse = TRUE, B = 6, na.rm = FALSE, seed = 5, cv_nfold = 3, gradient_params = list(eta = 0.3, max_depth = 2, subsample = 1, nrounds = 25),
    predict_sampled = "crossfit", predict_unsampled = "foldmean", bootstrap_refit = "leaves_only", boot_synthetic = synth, ...)))

test_that("boot_synthetic = TRUE leaves estimates, MSE and intervals bit-identical; the synthetic block has the documented shape", {
  skip_if_not_installed("xgboost"); M <- .lc_data()
  f0 <- .lc_fit(M, synth = FALSE); f1 <- .lc_fit(M)
  expect_identical(f1$ind, f0$ind); expect_identical(f1$var, f0$var); expect_identical(f1$CI, f0$CI)
  expect_null(f0$boot_diag$synthetic)
  s <- f1$boot_diag$synthetic; D <- nrow(f1$ind)
  expect_equal(length(s$domains), D); expect_equal(dim(s$boot), dim(s$truth)); expect_equal(nrow(s$boot), D); expect_equal(sum(s$sampled), 18)
  expect_equal(s$mse_synth, rowMeans((s$boot - s$truth)^2))
})
test_that("the replicate synthetic predictions equal a hand calculation (identity transformation, equal population weights)", {
  skip_if_not_installed("xgboost"); M <- .lc_data(); env <- new.env(); op <- options(povmap.leaves_only.diag_env = env); on.exit(options(op))
  f <- .lc_fit(M); s <- f$boot_diag$synthetic
  # hand calculation: with the identity back-transform and equal population weights the replicate's synthetic area prediction IS the replicate's booster
  # domain mean (the diag matrix of the leaves-only bootstrap), the area effect not included
  hand <- env$gb_pop_d_mean_boot_mat[match(s$domains, env$domains), , drop = FALSE]
  expect_equal(unname(s$boot), unname(hand), tolerance = 1e-12)
  # and the point fit's synthetic prediction of an unsampled area is its estimate (u = 0)
  un <- !s$sampled; expect_true(any(un)); expect_equal(s$point[un], f$ind$Mean[match(s$domains[un], as.character(f$ind$Domain))], tolerance = 1e-12)
})
test_that("the synthetic point prediction of an unsampled area equals its estimate also under the log transformation (back-transform with sigma_b^2)", {
  skip_if_not_installed("xgboost"); M <- .lc_data(); f <- .lc_fit(M, tr = "log", pw = "w"); s <- f$boot_diag$synthetic; un <- !s$sampled
  expect_equal(s$point[un], f$ind$Mean[match(s$domains[un], as.character(f$ind$Domain))], tolerance = 1e-10)
  expect_true(all(is.finite(s$boot[!is.na(s$boot)])))
  expect_error(suppressMessages(suppressWarnings(povmap::megb(fixed = y ~ x1 + x2, smp_data = M$smp, pop_data = M$pop, domains = "dom", transformation = "no", mse = TRUE, B = 2, na.rm = FALSE,
    bootstrap_refit = "lmm_only", predict_sampled = "full", predict_unsampled = "full", boot_synthetic = TRUE))), "leaves_only")
  expect_error(suppressMessages(povmap::megb(fixed = y ~ x1 + x2, smp_data = M$smp, pop_data = M$pop, domains = "dom", transformation = "no", mse = FALSE, na.rm = FALSE, boot_synthetic = TRUE)), "mse = TRUE")
})

## ---- megb_level_calibration on a hand-made fit ---------------------------------------------------------------------
.fake_fit <- function() {
  dom <- c("a", "b", "c", "d", "e")                      # a, b, c sampled; d, e unsampled
  list(boot_diag = list(synthetic = list(domains = dom, sampled = c(TRUE, TRUE, TRUE, FALSE, FALSE), point = c(0.50, 0.40, 0.30, 0.45, 0.35),
                                         boot = matrix(0, 5, 2), truth = matrix(0, 5, 2), mse_synth = c(0.010, 0.020, 0.030, 0.015, 0.025))),
       ind = data.frame(Domain = dom, Mean = c(.5, .4, .3, .45, .35), Mean_bench = c(.51, .41, .31, .46, .36)),
       var = data.frame(Domain = dom, Mean = c(.004, .005, .006, .007, .008), Var_bench = c(.0041, .0051, .0061, .0071, .0081)),
       CI = data.frame(Domain = dom, Lower = c(.4, .3, .2, .3, .2), Upper = c(.6, .5, .4, .6, .5), Lower_bench = c(.41, .31, .21, .31, .21), Upper_bench = c(.61, .51, .41, .61, .51)))
}
test_that("kappa_hat = M_hat / B_hat with the documented M_hat and B_hat; unsampled MSE scaled, sampled areas and estimates untouched", {
  f <- .fake_fit(); dir <- data.frame(Domain = c("a", "b", "c", "d"), ydir = c(0.60, 0.30, 0.50, 0.9), vdir = c(0.001, 0.002, 0.003, 0.001))   # d is unsampled: ignored
  r <- megb_level_calibration(f, dir, which = "bench")
  e2 <- c((0.60 - 0.50)^2 - 0.001, (0.30 - 0.40)^2 - 0.002, (0.50 - 0.30)^2 - 0.003); M <- mean(e2); B <- mean(c(0.010, 0.020, 0.030))
  expect_equal(r$M_hat, M); expect_equal(r$B_hat, B); expect_equal(r$kappa_raw, M / B); expect_equal(r$kappa, M / B); expect_equal(r$clipped, "none"); expect_false(r$fallback)
  expect_equal(r$n_areas, 3); tb <- r$table
  expect_identical(tb$MSE_D[1:3], tb$MSE_legacy[1:3]); expect_identical(tb$Lower_D[1:3], tb$Lower_legacy[1:3]); expect_identical(tb$Upper_D[1:3], tb$Upper_legacy[1:3])   # sampled areas untouched
  expect_equal(tb$MSE_D[4:5], (M / B) * c(.0071, .0081)); expect_equal(tb$Mean, f$ind$Mean_bench)
  z <- qnorm(0.975); expect_equal(tb$Upper_D[4], .46 + z * sqrt((M / B) * .0071))
})
test_that("kappa_hat is clipped to the bounds and falls back to 1 when M_hat or B_hat is not positive; precision weights change M_hat and B_hat", {
  f <- .fake_fit()
  big <- data.frame(Domain = c("a", "b", "c"), ydir = c(0.9, 0.0, 0.9), vdir = 0.001); r <- megb_level_calibration(f, big, which = "unbench")
  expect_equal(r$clipped, "upper"); expect_equal(r$kappa, 10); expect_gt(r$kappa_raw, 10)
  small <- data.frame(Domain = c("a", "b", "c"), ydir = c(0.5, 0.4, 0.3) + 1e-4, vdir = 0.001); r2 <- megb_level_calibration(f, small, which = "unbench")   # M_hat < 0: fallback
  expect_true(r2$fallback); expect_equal(r2$kappa, 1); expect_equal(r2$table$MSE_D, r2$table$MSE_legacy)
  low <- data.frame(Domain = c("a", "b", "c"), ydir = c(0.5, 0.4, 0.3) + 0.012, vdir = 0.0001); r3 <- megb_level_calibration(f, low, which = "unbench")
  expect_equal(r3$clipped, "lower"); expect_equal(r3$kappa, 0.1)
  dir <- data.frame(Domain = c("a", "b", "c"), ydir = c(0.60, 0.30, 0.50), vdir = c(0.001, 0.002, 0.003)); rp <- megb_level_calibration(f, dir, weights = "precision", which = "unbench")
  e2 <- c((0.60 - 0.50)^2 - 0.001, (0.30 - 0.40)^2 - 0.002, (0.50 - 0.30)^2 - 0.003); Bi <- c(0.010, 0.020, 0.030); a <- 1 / (c(0.001, 0.002, 0.003) + Bi)
  expect_equal(rp$M_hat, sum(a * e2) / sum(a)); expect_equal(rp$B_hat, sum(a * Bi) / sum(a))
  expect_error(megb_level_calibration(f, data.frame(Domain = "zzz", ydir = 0.5, vdir = 0.001)), "not in the fit")
})
test_that("megb_level_calibration works on a real fit; sampled areas unchanged; mse_type = 'var' rescales the half-widths by sqrt(kappa)", {
  skip_if_not_installed("xgboost"); M <- .lc_data(); f <- .lc_fit(M); s <- f$boot_diag$synthetic
  set.seed(1); dir <- data.frame(Domain = s$domains[s$sampled], ydir = s$point[s$sampled] + rnorm(sum(s$sampled), 0, 0.5), vdir = 0.01)
  r <- megb_level_calibration(f, dir, which = "unbench", mse_type = "mse"); tb <- r$table
  expect_true(is.finite(r$kappa)); expect_equal(tb$MSE_D[tb$sampled], tb$MSE_legacy[tb$sampled]); expect_equal(tb$MSE_D[!tb$sampled], r$kappa * tb$MSE_legacy[!tb$sampled])
  rv <- megb_level_calibration(f, dir, which = "unbench", mse_type = "var"); u <- !rv$table$sampled
  expect_equal(rv$table$Upper_D[u] - rv$table$Mean[u], sqrt(rv$kappa) * (rv$table$Upper_legacy[u] - rv$table$Mean[u]), tolerance = 1e-10)
})
