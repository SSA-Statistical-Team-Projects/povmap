# boot_synthetic and megb_level_calibration() (MEGB_BOOT_FIX_PREREG_v8_final.md, arm D, option A).
.lc_data <- function(D = 24, n_pop = 20, n_smp = 6, seed = 9, prop = FALSE) {
  set.seed(seed)
  pop <- data.frame(dom = factor(sprintf("d%02d", rep(seq_len(D), each = n_pop))), x1 = rnorm(D * n_pop), x2 = runif(D * n_pop), w = rpois(D * n_pop, 20) + 1)
  lp <- 0.3 * pop$x1 + rnorm(D, 0, 0.3)[as.integer(pop$dom)] + rnorm(nrow(pop), 0, 0.5)
  pop$y <- if (prop) plogis(lp) else exp(lp)
  keep <- unique(as.character(pop$dom))[1:18]
  smp <- do.call(rbind, lapply(keep, function(d) { z <- pop[pop$dom == d, ]; z[sample(nrow(z), n_smp), ] }))
  smp$wt <- runif(nrow(smp), 0.5, 2); rownames(smp) <- NULL
  list(pop = pop[, c("dom", "x1", "x2", "w")], smp = smp)
}
.lc_fit <- function(M, ..., tr = "no", pw = NULL, synth = TRUE, ps = "crossfit", pu = "foldmean", mse_type = "var")
  suppressMessages(suppressWarnings(povmap::megb(fixed = y ~ x1 + x2, smp_data = M$smp, smp_weights = "wt", pop_data = M$pop, pop_weights = pw, domains = "dom",
    transformation = tr, mse = TRUE, B = 6, na.rm = FALSE, seed = 5, cv_nfold = 3, gradient_params = list(eta = 0.3, max_depth = 2, subsample = 1, nrounds = 25),
    predict_sampled = ps, predict_unsampled = pu, bootstrap_refit = "leaves_only", mse_type = mse_type, boot_synthetic = synth, ...)))

test_that("boot_synthetic = TRUE leaves estimates, MSE and intervals bit-identical; the synthetic block has the documented shape", {
  skip_if_not_installed("xgboost"); M <- .lc_data()
  f0 <- .lc_fit(M, synth = FALSE); f1 <- .lc_fit(M)
  expect_identical(f1$ind, f0$ind); expect_identical(f1$var, f0$var); expect_identical(f1$CI, f0$CI)
  expect_null(f0$boot_diag$synthetic)
  s <- f1$boot_diag$synthetic; D <- nrow(f1$ind)
  expect_equal(length(s$domains), D); expect_equal(dim(s$boot), dim(s$truth)); expect_equal(nrow(s$boot), D); expect_equal(sum(s$sampled), 18)
  expect_equal(s$mse_synth, rowMeans((s$boot - s$truth)^2))
})

## ---- 1. the guard --------------------------------------------------------------------------------------------------
test_that("boot_synthetic = TRUE stops unless predict_sampled = 'crossfit'; other requirements stop too", {
  skip_if_not_installed("xgboost"); M <- .lc_data()
  expect_error(.lc_fit(M, ps = "full", pu = "foldmean"), "predict_sampled")
  expect_error(.lc_fit(M, ps = "full", pu = "full"), "predict_sampled")
  expect_error(suppressMessages(suppressWarnings(povmap::megb(fixed = y ~ x1 + x2, smp_data = M$smp, pop_data = M$pop, domains = "dom", transformation = "no", mse = TRUE, B = 2, na.rm = FALSE,
    bootstrap_refit = "lmm_only", predict_sampled = "full", predict_unsampled = "full", boot_synthetic = TRUE))), "predict_sampled|leaves_only")
  expect_error(suppressMessages(povmap::megb(fixed = y ~ x1 + x2, smp_data = M$smp, pop_data = M$pop, domains = "dom", transformation = "no", mse = FALSE, na.rm = FALSE, boot_synthetic = TRUE)), "mse = TRUE")
})

## ---- 2. held out in the replicates -----------------------------------------------------------------------------------
test_that("a sampled area's replicate synthetic prediction does not depend on its own replicate response (identity case)", {
  skip_if_not_installed("xgboost"); M <- .lc_data()
  base <- .lc_fit(M)$boot_diag$synthetic
  k <- "d03"; stopifnot(k %in% base$domains[base$sampled])
  op <- options(povmap.boot_synthetic.perturb = list(domain = k, delta = 7)); on.exit(options(op))
  pert <- .lc_fit(M)$boot_diag$synthetic; options(op)
  ik <- match(k, base$domains)
  expect_identical(pert$boot[ik, ], base$boot[ik, ])                                   # area k: unchanged by its own replicate data
  expect_gt(max(abs(pert$boot[-ik, ] - base$boot[-ik, ])), 1e-8)                        # the perturbation does reach other areas' predictions (the test can fail)
  expect_identical(pert$point, base$point)                                              # the point fit is untouched
})

## ---- 3. the transformed path: per-cell back-transform with the replicate's sigma_b*^2 and no area effect -----------------
.hand_syn <- function(f, env, M, bt) {
  s <- f$boot_diag$synthetic; dom <- as.character(M$pop$dom); w <- M$pop$w; s2e <- f$boot_diag$error_sd_orig^2
  vapply(seq_along(env$gb_pop_cell_boot_list), function(b) {
    h <- bt(env$gb_pop_cell_boot_list[[b]], s2e + env$sb2_boot_vec[b])
    a <- tapply(h * w, dom, sum) / tapply(w, dom, sum); as.numeric(a[match(s$domains, names(a))]) }, numeric(length(s$domains)))
}
test_that("log: boot = per-cell exp(x + (sigma_e^2 + sigma_b*^2) / 2), aggregated with the population weights, no area effect", {
  skip_if_not_installed("xgboost"); M <- .lc_data(); env <- new.env(); op <- options(povmap.leaves_only.diag_env = env); on.exit(options(op))
  f <- .lc_fit(M, tr = "log", pw = "w"); s <- f$boot_diag$synthetic
  hand <- .hand_syn(f, env, M, function(x, v) exp(x + 0.5 * v))
  expect_equal(unname(s$boot), unname(hand), tolerance = 1e-10)
})
test_that("arcsin (plain and exact): boot = per-cell sin^2 / exact E[sin^2] with the replicate's sigma_b*^2, no area effect", {
  skip_if_not_installed("xgboost"); M <- .lc_data(prop = TRUE); env <- new.env(); op <- options(povmap.leaves_only.diag_env = env); on.exit(options(op))
  f <- .lc_fit(M, tr = "arcsin", pw = "w", arcsin_backtransform = "plain"); s <- f$boot_diag$synthetic
  expect_equal(unname(s$boot), unname(.hand_syn(f, env, M, function(x, v) sin(pmax(0, pmin(x, pi / 2)))^2)), tolerance = 1e-10)
  env2 <- new.env(); options(povmap.leaves_only.diag_env = env2)
  fe <- .lc_fit(M, tr = "arcsin", pw = "w", arcsin_backtransform = "exact"); se <- fe$boot_diag$synthetic
  expect_equal(unname(se$boot), unname(.hand_syn(fe, env2, M, function(x, v) { m <- pmax(0, pmin(x, pi / 2)); 0.5 - 0.5 * cos(2 * m) * exp(-2 * v) })), tolerance = 1e-10)
  expect_gt(max(abs(se$boot - s$boot)), 1e-6)                                            # the exact form differs from the plain one: the variance enters
})
test_that("identity case: boot equals the replicate booster domain mean (hand calculation), and the point synthetic of an unsampled area is its estimate", {
  skip_if_not_installed("xgboost"); M <- .lc_data(); env <- new.env(); op <- options(povmap.leaves_only.diag_env = env); on.exit(options(op))
  f <- .lc_fit(M); s <- f$boot_diag$synthetic
  expect_equal(unname(s$boot), unname(env$gb_pop_d_mean_boot_mat[match(s$domains, env$domains), , drop = FALSE]), tolerance = 1e-12)
  un <- !s$sampled; expect_true(any(un)); expect_equal(s$point[un], f$ind$Mean[match(s$domains[un], as.character(f$ind$Domain))], tolerance = 1e-12)
  fl <- .lc_fit(M, tr = "log", pw = "w"); sl <- fl$boot_diag$synthetic; u2 <- !sl$sampled
  expect_equal(sl$point[u2], fl$ind$Mean[match(sl$domains[u2], as.character(fl$ind$Domain))], tolerance = 1e-10)
})

## ---- 4. held-out against the full fit --------------------------------------------------------------------------------
test_that("a sampled area's stored point synthetic prediction is the held-out fold model's prediction and differs from the full-fit booster's", {
  skip_if_not_installed("xgboost"); M <- .lc_data(); f <- .lc_fit(M); s <- f$boot_diag$synthetic
  ff <- povmap:::.megb_fold_fit(f$megb_model); skip_if(is.null(ff), "no fold models in the fit")
  X <- data.matrix(M$pop[, c("x1", "x2")]); dom <- as.character(M$pop$dom)
  slice <- function(b) { bi <- xgboost::xgb.attr(b, "best_iteration"); nr <- if (is.null(bi)) xgboost::xgb.get.num.boosted.rounds(b) else as.integer(bi) + 1L; xgboost::xgb.slice.Booster(b, 1L, nr) }
  held <- full <- rep(NA_real_, length(s$domains))
  for (k in seq_along(ff$fold_models)) for (d in ff$fold_domains[[k]]) {
    i <- match(d, s$domains); r <- dom == d
    held[i] <- mean(predict(slice(ff$fold_models[[k]]), X[r, , drop = FALSE])); full[i] <- mean(predict(f$model, X[r, , drop = FALSE]))
  }
  sm <- which(s$sampled); expect_false(anyNA(held[sm]))
  expect_equal(s$point[sm], held[sm], tolerance = 1e-6)                                  # the held-out fold model's prediction (identity transformation, equal pop weights)
  expect_gt(min(abs(s$point[sm] - full[sm])), 0); expect_gt(mean(abs(s$point[sm] - full[sm]) > 1e-8), 0.9)   # and it is not the full booster's
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

## ---- 6. n_left_out ---------------------------------------------------------------------------------------------------
test_that("n_left_out counts every sampled area left out of M_hat and B_hat, also those missing from direct", {
  f <- .fake_fit()
  d2 <- data.frame(Domain = c("a", "b"), ydir = c(0.6, 0.3), vdir = c(0.001, 0.002)); r <- megb_level_calibration(f, d2, which = "unbench")
  expect_equal(r$n_areas, 2); expect_equal(r$n_left_out, 1)                              # c is sampled but absent from direct
  d3 <- data.frame(Domain = c("a", "b", "c"), ydir = c(0.6, NA, 0.5), vdir = c(0.001, 0.002, NA)); r3 <- megb_level_calibration(f, d3, which = "unbench")
  expect_equal(r3$n_areas, 1); expect_equal(r3$n_left_out, 2)                            # b: NA ydir; c: NA vdir
  d4 <- data.frame(Domain = c("a", "b", "c", "d"), ydir = c(0.6, 0.3, 0.5, 0.9), vdir = 0.001); r4 <- megb_level_calibration(f, d4, which = "unbench")
  expect_equal(r4$n_left_out, 0)                                                         # d is unsampled: not counted
})

## ---- 5. kappa = 1 gives the legacy intervals ---------------------------------------------------------------------------
test_that("with kappa_hat = 1 (fallback) the rebuilt unsampled intervals equal the legacy ones, for mse_type 'mse' (same clamp) and 'var'", {
  skip_if_not_installed("xgboost")
  force1 <- function(s) data.frame(Domain = s$domains[s$sampled], ydir = s$point[s$sampled], vdir = 0.01)     # ydir = synthetic point: M_hat < 0, so the fallback sets kappa_hat = 1
  M <- .lc_data(); fv <- .lc_fit(M, mse_type = "var"); rv <- megb_level_calibration(fv, force1(fv$boot_diag$synthetic), which = "unbench", mse_type = "var")
  expect_true(rv$fallback); expect_equal(rv$kappa, 1)
  tb <- rv$table; expect_identical(tb$Lower_D, fv$CI$Lower[match(tb$Domain, fv$CI$Domain)]); expect_identical(tb$Upper_D, fv$CI$Upper[match(tb$Domain, fv$CI$Domain)])
  fm <- .lc_fit(M, mse_type = "mse"); rm <- megb_level_calibration(fm, force1(fm$boot_diag$synthetic), which = "unbench", mse_type = "mse"); tm <- rm$table
  expect_true(rm$fallback); expect_identical(tm$Lower_D, fm$CI$Lower[match(tm$Domain, fm$CI$Domain)]); expect_identical(tm$Upper_D, fm$CI$Upper[match(tm$Domain, fm$CI$Domain)])
  Mp <- .lc_data(prop = TRUE); fa <- .lc_fit(Mp, tr = "arcsin", pw = "w", mse_type = "mse", arcsin_backtransform = "plain")        # bounded outcome: the fit clamps the interval to [0, 1]
  ra <- megb_level_calibration(fa, force1(fa$boot_diag$synthetic), which = "unbench", mse_type = "mse", clamp = c(0, 1)); ta <- ra$table
  expect_true(ra$fallback); expect_identical(ta$Lower_D, fa$CI$Lower[match(ta$Domain, fa$CI$Domain)]); expect_identical(ta$Upper_D, fa$CI$Upper[match(ta$Domain, fa$CI$Domain)])
})
test_that("megb_level_calibration works on a real fit; sampled areas unchanged; mse_type = 'var' rescales the half-widths by sqrt(kappa)", {
  skip_if_not_installed("xgboost"); M <- .lc_data(); f <- .lc_fit(M); s <- f$boot_diag$synthetic
  set.seed(1); dir <- data.frame(Domain = s$domains[s$sampled], ydir = s$point[s$sampled] + rnorm(sum(s$sampled), 0, 0.5), vdir = 0.01)
  r <- megb_level_calibration(f, dir, which = "unbench", mse_type = "mse"); tb <- r$table
  expect_true(is.finite(r$kappa)); expect_equal(tb$MSE_D[tb$sampled], tb$MSE_legacy[tb$sampled]); expect_equal(tb$MSE_D[!tb$sampled], r$kappa * tb$MSE_legacy[!tb$sampled])
  rv <- megb_level_calibration(f, dir, which = "unbench", mse_type = "var"); u <- !rv$table$sampled
  expect_equal(rv$table$Upper_D[u] - rv$table$Mean[u], sqrt(rv$kappa) * (rv$table$Upper_legacy[u] - rv$table$Mean[u]), tolerance = 1e-10)
})
