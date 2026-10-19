# boot_area_var (MEGB_BOOT_FIX_PREREG_v3.md): "legacy" is the previous behaviour bit for bit; "centred" and "centred_mm" change only the bootstrap world.
.bav_data <- function(D = 30, n_pop = 20, n_smp = 8, seed = 7, single = TRUE) {
  set.seed(seed)
  pop <- data.frame(dom = factor(rep(sprintf("d%02d", 1:D), each = n_pop)),
                    x1 = rnorm(D * n_pop), x2 = runif(D * n_pop), w = rpois(D * n_pop, 20) + 1)
  pop$xm <- rnorm(D)[as.integer(pop$dom)]
  pop$y <- 0.5 * pop$x1 + 0.4 * pop$xm + rnorm(D, 0, 0.3)[as.integer(pop$dom)] + rnorm(nrow(pop))
  doms <- split(pop, pop$dom)[1:25]
  smp <- do.call(rbind, lapply(seq_along(doms), function(i) { d <- doms[[i]]; d[sample(nrow(d), if (single && i %in% c(3, 11)) 1L else n_smp), ] }))
  smp$wt <- runif(nrow(smp), 0.3, 3)
  list(pop = pop, smp = smp)
}
.bav_fit <- function(M, mode, weightedBS = FALSE, bootstrap_refit = "leaves_only", lib = NULL, B = 4, ps = "crossfit", pu = "foldmean")
  suppressMessages(suppressWarnings(povmap::megb(
    fixed = y ~ x1 + x2 + xm, smp_data = M$smp, smp_weights = "wt",
    pop_data = M$pop[, c("dom", "x1", "x2", "xm", "w")], pop_weights = "w", domains = "dom",
    transformation = "no", mse = TRUE, B = B, na.rm = FALSE, seed = 5, cv_nfold = 5, weightedBS = weightedBS,
    gradient_params = list(eta = 0.3, max_depth = 3, subsample = 1, nrounds = 60),
    predict_sampled = ps, predict_unsampled = pu, bootstrap_refit = bootstrap_refit,
    boot_area_var = mode)))

## ---- toy dataset with known answers -------------------------------------------------------
## three sampled areas: a (3 units), b (2 units), c (1 unit); unit weights 1
.bav_toy <- function() {
  smp <- data.frame(dom = c("a", "a", "a", "b", "b", "c"))
  list(smp = smp, Y = c(1, 2, 6, 10, 14, 20), gbs = c(0, 0, 0, 8, 8, 17), res = c(1, 2, 6, 2, -2, 3), w = rep(1, 6))
}

test_that("centred: residuals are centred within the area, scaled by sqrt(n/(n-1)); the single-unit area is excluded from the pool", {
  T0 <- .bav_toy()
  out <- povmap:::.megb_boot_area_var("centred", T0$res, T0$Y, T0$gbs, T0$smp, "dom", error_sd = 1, weights = T0$w, weightedBS = FALSE, ran_effs_legacy = c(.1, -.1, 0))
  ea <- c(1, 2, 6) - 3                        # area a: centred (mean 3), scaled by sqrt(3/2)
  eb <- c(2, -2) - 0                          # area b: centred (mean 0), scaled by sqrt(2/1)
  expect_equal(out$gb_res[1:3], ea * sqrt(3 / 2)); expect_equal(out$gb_res[4:5], eb * sqrt(2))
  expect_equal(out$pool_res, c(ea * sqrt(3 / 2), eb * sqrt(2)))      # 5 values: area c (one unit) is not in the pool
  expect_null(out$pool_w); expect_equal(out$ran_effs, c(.1, -.1, 0)); expect_identical(names(out$area_mm), c("pool_sd", "error_sd"))    # u* as legacy; only the pool diagnostics
})

test_that("centred with weightedBS uses the weighted within-area mean and weights the pool", {
  T0 <- .bav_toy(); w <- c(1, 1, 2, 1, 3, 1)
  out <- povmap:::.megb_boot_area_var("centred", T0$res, T0$Y, T0$gbs, T0$smp, "dom", 1, w, TRUE, c(0, 0, 0))
  mu_a <- sum(c(1, 2, 6) * c(1, 1, 2)) / 4; mu_b <- sum(c(2, -2) * c(1, 3)) / 4
  expect_equal(out$gb_res[1:3], (c(1, 2, 6) - mu_a) * sqrt(3 / 2)); expect_equal(out$gb_res[4:5], (c(2, -2) - mu_b) * sqrt(2))
  expect_equal(out$pool_w, c(1, 1, 2, 1, 3))
})

test_that("centred_mm: sigma_b2(MM), the centred deviations d_i and their rescale are the documented moment estimator", {
  T0 <- .bav_toy(); sd_e <- 0.3            # small error_sd so that sigma_b2(MM) > 0 and the checks below run
  out <- povmap:::.megb_boot_area_var("centred_mm", T0$res, T0$Y, T0$gbs, T0$smp, "dom", sd_e, T0$w, FALSE, c(0, 0, 0))
  d <- c(mean(c(1, 2, 6)), mean(c(2, 6)), 3)           # Y - gb_smp by area: a (1,2,6), b (2,6), c (3); the unit weights are 1
  a <- c(3, 2, 1); n <- a; dbar <- sum(a * d) / sum(a); vw <- sum(a * (d - dbar)^2) / sum(a); nz <- sum(a * sd_e^2 / n) / sum(a)
  sb2 <- max(0, vw - nz)
  expect_equal(out$area_mm$sigma_b2_mm, sb2); expect_equal(out$area_mm$var_w_d, vw); expect_equal(out$area_mm$mean_noise, nz)
  expect_equal(out$area_mm$n_single_unit_areas, 1)
  expect_gt(sb2, 0); expect_gt(out$area_mm$sigma_b2_mm, 0)
  expect_equal(sum(a * out$ran_effs) / sum(a), 0, tolerance = 1e-12)                       # centred with the area weights
  expect_equal(sqrt(sum(a * out$ran_effs^2) / sum(a)), sqrt(sb2), tolerance = 1e-12)       # rescaled to sd sigma_b(MM)
  expect_equal(out$ran_effs, (d - dbar) / sqrt(vw) * sqrt(sb2))
})

test_that("centred_mm: sigma_b2(MM) is clipped at zero and then there is no area effect", {
  T0 <- .bav_toy()
  out <- povmap:::.megb_boot_area_var("centred_mm", T0$res, T0$Y, T0$gbs, T0$smp, "dom", error_sd = 50, weights = T0$w, weightedBS = FALSE, ran_effs_legacy = c(1, 2, 3))
  expect_equal(out$area_mm$sigma_b2_mm, 0); expect_equal(out$ran_effs, c(0, 0, 0))
})

test_that("block_sample: single-unit and unsampled areas draw from the pool; legacy arguments change nothing", {
  smp <- data.frame(dom = c("a", "a", "b", "c")); gb <- c(1, 2, 3, 99); domains <- c("a", "b", "c", "d"); in_samp <- c(TRUE, TRUE, TRUE, FALSE)
  pop <- data.frame(dom = c("a", "a", "b", "c", "d", "d", "d"))
  set.seed(1); leg <- povmap:::block_sample(domains, in_samp, smp, "dom", smp, gb)
  set.seed(1); leg2 <- povmap:::block_sample(domains, in_samp, smp, "dom", smp, gb, pool = NULL, pool_w = NULL); expect_identical(leg, leg2)
  set.seed(2); x <- povmap:::block_sample(domains, in_samp, smp, "dom", pop, gb, pool = c(10, 11))
  expect_true(all(x[1:2] %in% gb[1:2]))                 # area a keeps its own residuals
  expect_true(all(x[3:7] %in% c(10, 11)))               # b (one unit), c (one unit) and the unsampled d draw from the pool
})

## ---- integration -----------------------------------------------------------------------------
test_that("the three options give the same point estimates and legacy is the default", {
  skip_if_not_installed("xgboost"); M <- .bav_data()
  f0 <- .bav_fit(M, "legacy"); f1 <- .bav_fit(M, "centred"); f2 <- .bav_fit(M, "centred_mm")
  fd <- suppressMessages(suppressWarnings(povmap::megb(fixed = y ~ x1 + x2 + xm, smp_data = M$smp, smp_weights = "wt", pop_data = M$pop[, c("dom", "x1", "x2", "xm", "w")], pop_weights = "w",
    domains = "dom", transformation = "no", mse = TRUE, B = 4, na.rm = FALSE, seed = 5, cv_nfold = 5, weightedBS = FALSE, gradient_params = list(eta = 0.3, max_depth = 3, subsample = 1, nrounds = 60),
    predict_sampled = "crossfit", predict_unsampled = "foldmean")))
  expect_identical(f1$ind$Mean, f0$ind$Mean); expect_identical(f2$ind$Mean, f0$ind$Mean)
  expect_identical(fd$ind, f0$ind); expect_identical(fd$var, f0$var); expect_identical(fd$CI, f0$CI)       # default = legacy, bit for bit
  expect_false(isTRUE(all.equal(f1$var$Mean, f0$var$Mean))); expect_false(isTRUE(all.equal(f2$var$Mean, f1$var$Mean)))   # the bootstrap does change
  expect_null(f0$boot_diag$boot_area_var); expect_equal(f2$boot_diag$boot_area_var$mode, "centred_mm"); expect_true(is.finite(f2$boot_diag$boot_area_var$sigma_b2_mm))
  expect_equal(f1$boot_diag$boot_area_var$n_pool, sum(table(M$smp$dom) >= 2) * 0 + sum(as.numeric(table(M$smp$dom))[as.numeric(table(M$smp$dom)) >= 2]))
})

test_that("the options also run with weightedBS and the full / lmm_only refits", {
  skip_if_not_installed("xgboost"); M <- .bav_data()
  for (m in c("centred", "centred_mm")) for (r in c("full", "lmm_only")) {
    lm <- identical(r, "lmm_only")     # lmm_only has no fold-model bootstrap, so "crossfit"/"foldmean" are unavailable there
    if (lm && m == "centred_mm") {     # lmm_only keeps no fold models, so there is no held-out prediction for the d_i: an error, not an in-sample fallback
      expect_error(.bav_fit(M, m, weightedBS = TRUE, bootstrap_refit = r, B = 2, ps = "full", pu = "full"), "cross-fitted")
      next
    }
    f <- .bav_fit(M, m, weightedBS = TRUE, bootstrap_refit = r, B = 2, ps = if (lm) "full" else "crossfit", pu = if (lm) "full" else "foldmean"); expect_true(all(is.finite(f$var$Mean[!is.na(f$var$Mean)])))
  }
})

test_that("legacy is bit-identical to povmap 1.0.7 (tag xgb-convention-v0.6.5)", {
  skip_if_not_installed("xgboost"); lib <- Sys.getenv("POVMAP_LIB_107", "D:/David/R/lib-povmap-9888899"); skip_if_not(dir.exists(file.path(lib, "povmap")), "povmap 1.0.7 library not available")
  M <- .bav_data(); here <- .bav_fit(M, "legacy", B = 3)
  tf <- tempfile(fileext = ".rds"); sf <- tempfile(fileext = ".R"); saveRDS(M, tf)
  writeLines(c(sprintf(".libPaths(c('%s', .libPaths()))", lib), "library(povmap); M <- readRDS(commandArgs(TRUE)[1])",
    "f <- suppressMessages(suppressWarnings(megb(fixed = y ~ x1 + x2 + xm, smp_data = M$smp, smp_weights = 'wt', pop_data = M$pop[, c('dom','x1','x2','xm','w')], pop_weights = 'w', domains = 'dom', transformation = 'no', mse = TRUE, B = 3, na.rm = FALSE, seed = 5, cv_nfold = 5, gradient_params = list(eta = 0.3, max_depth = 3, subsample = 1, nrounds = 60), predict_sampled = 'crossfit', predict_unsampled = 'foldmean', bootstrap_refit = 'leaves_only', weightedBS = FALSE)))",
    "saveRDS(list(ind = f$ind, var = f$var, CI = f$CI, bd = f$boot_diag, v = as.character(packageVersion('povmap'))), commandArgs(TRUE)[2])"), sf)
  of <- tempfile(fileext = ".rds"); rs <- file.path(R.home("bin"), "Rscript"); system2(rs, c(sf, tf, of), stdout = FALSE, stderr = FALSE)
  skip_if_not(file.exists(of), "could not run the 1.0.7 reference"); old <- readRDS(of); expect_equal(old$v, "1.0.7")
  expect_identical(here$ind, old$ind); expect_identical(here$var, old$var); expect_identical(here$CI, old$CI)
  expect_identical(here$boot_diag$boot_ran_eff_sd, old$bd$boot_ran_eff_sd); expect_identical(here$boot_diag$boot_error_sd, old$bd$boot_error_sd)
})

## ---- v4 changes (MEGB_BOOT_FIX_PREREG_v4.md): Kish noise term, drawn-u* sd, diagnostics ------------------------
test_that("centred_mm: with unit weights the noise term uses the Kish size (sum w)^2 / sum w^2; without weights n_i", {
  T0 <- .bav_toy(); sd_e <- 0.3; w <- c(1, 3, 1, 2, 2, 5)           # areas a (1,3,1), b (2,2), c (5)
  a <- c(5, 4, 5); n <- c(3, 2, 1); neff <- c(25 / 11, 16 / 8, 25 / 25)
  d <- c((1 * 0 + 3 * 0 + 1 * 0) + c(1 * 1 + 3 * 2 + 1 * 6) / 5 - 0, (2 * 2 + 2 * 6) / 4, 3)       # weighted mean of Y - gb_smp by area
  dbar <- sum(a * d) / sum(a); vw <- sum(a * (d - dbar)^2) / sum(a)
  out <- povmap:::.megb_boot_area_var("centred_mm", T0$res, T0$Y, T0$gbs, T0$smp, "dom", sd_e, w, FALSE, c(0, 0, 0))
  expect_equal(out$area_mm$mean_noise, sum(a * sd_e^2 / neff) / sum(a)); expect_equal(out$area_mm$mean_noise_n, sum(a * sd_e^2 / n) / sum(a))
  expect_equal(out$area_mm$var_w_d, vw); expect_equal(out$area_mm$sigma_b2_mm, max(0, vw - sum(a * sd_e^2 / neff) / sum(a)))
  expect_equal(out$area_mm$error_sd, sd_e); expect_true(is.finite(out$area_mm$pool_sd))
  ## the draws are equal-probability when weightedBS is FALSE: the pool is then centred and rescaled with equal weights
  expect_gt(out$area_mm$sigma_b2_mm, 0); expect_equal(mean(out$ran_effs), 0, tolerance = 1e-12); expect_equal(sqrt(mean(out$ran_effs^2)), sqrt(out$area_mm$sigma_b2_mm), tolerance = 1e-12)
  o2 <- povmap:::.megb_boot_area_var("centred_mm", T0$res, T0$Y, T0$gbs, T0$smp, "dom", sd_e, w, TRUE, c(0, 0, 0))     # weightedBS: draws with probability a_i
  expect_gt(o2$area_mm$sigma_b2_mm, 0); expect_equal(sum(a * o2$ran_effs) / sum(a), 0, tolerance = 1e-12); expect_equal(sqrt(sum(a * o2$ran_effs^2) / sum(a)), sqrt(o2$area_mm$sigma_b2_mm), tolerance = 1e-12)
  o3 <- povmap:::.megb_boot_area_var("centred_mm", T0$res, T0$Y, T0$gbs, T0$smp, "dom", sd_e, NULL, FALSE, c(0, 0, 0))   # no weights: n_i
  expect_equal(o3$area_mm$mean_noise, o3$area_mm$mean_noise_n)
})

.bav_data2 <- function() { M <- .bav_data(D = 40, n_smp = 10, seed = 11, single = FALSE); set.seed(3); a <- rnorm(40, 0, 0.7); M$pop$y <- M$pop$y + a[as.integer(M$pop$dom)]
  M$smp <- do.call(rbind, lapply(split(M$pop, M$pop$dom)[1:30], function(d) d[sample(nrow(d), 10), ])); M$smp$wt <- runif(nrow(M$smp), 0.3, 3); M }
test_that("the empirical sd of the drawn u* equals sqrt(sigma_b2(MM)) under centred_mm and ran_eff_sd under legacy", {
  skip_if_not_installed("xgboost"); M <- .bav_data2()
  draw_sd <- function(mode, wbs) { env <- new.env(); op <- options(povmap.boot_area_var.diag_env = env); on.exit(options(op))
    f <- .bav_fit(M, mode, weightedBS = wbs, B = 400); u <- env$u_d_star; list(f = f, sd = sd(as.vector(u)), m = mean(as.vector(u)), ncol = ncol(u)) }
  for (wbs in c(FALSE, TRUE)) {
    g <- draw_sd("centred_mm", wbs); s2 <- g$f$boot_diag$boot_area_var$sigma_b2_mm; expect_true(s2 > 0); expect_equal(g$ncol, 400)
    expect_equal(g$sd, sqrt(s2), tolerance = 0.05)                                    # Monte Carlo tolerance: 400 x 40 draws
  }
  l <- draw_sd("legacy", TRUE); expect_equal(l$sd, l$f$megb_model$ran_eff_sd, tolerance = 0.06)       # legacy: weighted rescale to ran_eff_sd, draws with the same probabilities
})

test_that("the diagnostics carry the pool sd, the error sd and the Kish-adjusted mean noise", {
  skip_if_not_installed("xgboost"); M <- .bav_data(); f <- .bav_fit(M, "centred_mm"); b <- f$boot_diag$boot_area_var
  expect_true(all(c("pool_sd", "error_sd", "mean_noise", "mean_noise_n", "sigma_b2_mm", "var_w_d", "n_pool", "mode") %in% names(b)))
  expect_equal(b$error_sd, f$boot_diag$error_sd_orig); expect_true(b$mean_noise >= b$mean_noise_n * (1 - 1e-12))     # the Kish size is never larger than n_i
  fc <- .bav_fit(M, "centred"); expect_true(all(c("pool_sd", "error_sd", "n_pool") %in% names(fc$boot_diag$boot_area_var)))
})

## ---- v4 implementation correction: d_i from held-out (cross-fitted) predictions, whatever predict_sampled is ----
test_that("centred_mm: d_i uses the held-out fold predictions, also with predict_sampled = 'full'", {
  skip_if_not_installed("xgboost"); M <- .bav_data2()
  run <- function(ps) { env <- new.env(); op <- options(povmap.boot_area_var.diag_env = env); on.exit(options(op))
    f <- .bav_fit(M, "centred_mm", B = 2, ps = ps, pu = "foldmean"); list(f = f, e = env) }
  sb2 <- function(r, g) { w <- M$smp$wt[order(as.character(M$smp$dom))]; fd <- factor(r$e$dom); a <- as.numeric(table(fd)); d <- as.numeric(tapply(w * (r$e$Y - g), fd, sum) / tapply(w, fd, sum))   # megb rescales the unit weights within each area to sum n_i
    neff <- as.numeric(tapply(w, fd, sum))^2 / as.numeric(tapply(w^2, fd, sum)); sdE <- r$f$boot_diag$error_sd_orig; dbar <- sum(a * d) / sum(a); max(0, sum(a * (d - dbar)^2) / sum(a) - sum(a * sdE^2 / neff) / sum(a)) }
  for (ps in c("crossfit", "full")) {
    r <- run(ps)
    expect_false(isTRUE(all.equal(r$e$gb_cf, r$e$gb_smp)) && ps == "full")                # 'full': gb_smp is the in-sample fit, gb_cf the held-out one
    expect_equal(r$f$boot_diag$boot_area_var$sigma_b2_mm, sb2(r, r$e$gb_cf), tolerance = 1e-8)   # sigma_b2(MM) is computed from the held-out prediction
    if (ps == "full") expect_false(isTRUE(all.equal(r$f$boot_diag$boot_area_var$sigma_b2_mm, sb2(r, r$e$gb_smp), tolerance = 1e-6)))   # not from the in-sample one
  }
})
