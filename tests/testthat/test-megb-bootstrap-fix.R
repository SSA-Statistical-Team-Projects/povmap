# megb's leaves_only bootstrap (3 Oct 2026):
#  * the leaf refresh carries the rescaled survey weights the booster was trained
#    with (it was unweighted);
#  * each replicate's random effects are re-estimated as the point estimate's are:
#    by ML, from the out-of-fold residuals of the fold models (each refreshed on
#    its own training rows), not by REML from the refreshed full booster's
#    in-sample residuals.
.bf_data <- function(D = 30, n_pop = 20, n_smp = 8, seed = 7) {
  set.seed(seed)
  pop <- data.frame(dom = factor(rep(sprintf("d%02d", 1:D), each = n_pop)),
                    x1 = rnorm(D * n_pop), x2 = runif(D * n_pop), w = rpois(D * n_pop, 20) + 1)
  pop$xm <- rnorm(D)[as.integer(pop$dom)]
  pop$y <- 0.5 * pop$x1 + 0.4 * pop$xm + rnorm(D, 0, 0.3)[as.integer(pop$dom)] + rnorm(nrow(pop))
  smp <- do.call(rbind, lapply(split(pop, pop$dom)[1:25], function(d) d[sample(nrow(d), n_smp), ]))
  smp$wt <- runif(nrow(smp), 0.3, 3)                # unequal survey weights
  list(pop = pop, smp = smp)
}
.bf_params <- list(eta = 0.3, max_depth = 3, subsample = 1, nrounds = 100)
## The fix is tested in the configuration the exhibits ran (five folds, full prediction) unless a test sets otherwise;
## the defaults from povmap 1.0.1 are tested in test-megb-defaults.R.
.bf_fit <- function(M, mse = TRUE, B = 5, cv_nfold = 5, predict_sampled = "full", predict_unsampled = "full", ...)
  suppressMessages(suppressWarnings(povmap::megb(
  fixed = y ~ x1 + x2 + xm, smp_data = M$smp, smp_weights = "wt",
  pop_data = M$pop[, c("dom", "x1", "x2", "xm", "w")], pop_weights = "w", domains = "dom",
  transformation = "no", mse = mse, B = B, na.rm = FALSE, seed = 5,
  gradient_params = .bf_params, cv_nfold = cv_nfold, predict_sampled = predict_sampled,
  predict_unsampled = predict_unsampled, ...)))
.rescale <- function(w, d) { s <- ave(w, d, FUN = sum); n <- ave(w, d, FUN = length); w * n / s }

test_that("the weighted leaf refresh reproduces a booster on its own data; the unweighted one does not", {
  skip_if_not_installed("xgboost")
  M <- .bf_data(); X <- M$smp[, c("x1", "x2", "xm")]; w <- .rescale(M$smp$wt, M$smp$dom)
  r <- povmap:::train_gbmodel("xgboost", X, M$smp$y, params = .bf_params, weights = w,
                              groups = M$smp$dom, early_stopping = FALSE)
  Xm <- data.matrix(X); n <- xgboost::xgb.get.num.boosted.rounds(r$boosting)
  p0 <- predict(r$boosting, Xm)
  expect_identical(predict(povmap:::.megb_refresh(r$boosting, Xm, M$smp$y, w, n), Xm), p0)
  expect_gt(max(abs(predict(povmap:::.megb_refresh(r$boosting, Xm, M$smp$y, NULL, n), Xm) - p0)), 1e-3)
})

test_that("on unperturbed data the bootstrap's random-effect step reproduces the point estimate's", {
  skip_if_not_installed("xgboost")
  M <- .bf_data()
  smp <- M$smp[, c("dom", "x1", "x2", "xm")]; pop <- M$pop[, c("dom", "x1", "x2", "xm")]
  f <- suppressMessages(povmap:::megb_em(Y = M$smp$y, X = smp[, -1], dom_name = "dom", smp_data = smp,
         pop_data = pop, gradient_params = .bf_params, seed = 5, smp_weights_vec = M$smp$wt,
         cv_folds = "domain", keep_fold_models = TRUE))
  mm <- f$megb_model
  ## the bootstrap's view: rows sorted by domain, weights rescaled within domain
  sd <- f$inp_smp_data$smp_data; Y <- f$inp_smp_data$target_var
  so <- order(as.character(sd$dom)); sds <- sd[so, , drop = FALSE]
  Xs <- as.matrix(sds[, f$cov_names_proc, drop = FALSE]); ws <- .rescale(M$smp$wt[so], sds$dom)
  CFB <- povmap:::.megb_fold_boot_setup(povmap:::.megb_fold_fit(mm), Xs, so, as.character(sds$dom),
                                        NULL, NULL, NULL, NULL, use_leaf = FALSE, predict = FALSE)
  st  <- povmap:::.megb_fold_boot_step(CFB, Xs, NULL, Y[so], ws, refresh = FALSE)
  expect_identical(st$oof, mm$oof_prediction[so])
  fit_data <- sds
  lm  <- povmap:::.megb_oof_re_fit(fit_data, r ~ 1 + (1 | dom), Y[so] - st$oof, ws)
  u0  <- lme4::ranef(mm$effect_model)$dom; u1 <- lme4::ranef(lm)$dom
  expect_equal(u1[rownames(u0), 1], u0[, 1], tolerance = 1e-6)
  expect_equal(as.data.frame(lme4::VarCorr(lm))$sdcor, as.data.frame(lme4::VarCorr(mm$effect_model))$sdcor, tolerance = 1e-6)
  expect_false(lme4::isREML(lm))
})

test_that("the fold setup refuses fold models that do not reproduce the fit's out-of-fold predictions", {
  skip_if_not_installed("xgboost")
  M <- .bf_data()
  smp <- M$smp[, c("dom", "x1", "x2", "xm")]; pop <- M$pop[, c("dom", "x1", "x2", "xm")]
  f <- suppressMessages(povmap:::megb_em(Y = M$smp$y, X = smp[, -1], dom_name = "dom", smp_data = smp,
         pop_data = pop, gradient_params = .bf_params, seed = 5, smp_weights_vec = M$smp$wt,
         cv_folds = "domain", keep_fold_models = TRUE))
  ff <- povmap:::.megb_fold_fit(f$megb_model); sd <- f$inp_smp_data$smp_data
  Xs <- as.matrix(sd[, f$cov_names_proc, drop = FALSE])
  so <- order(as.character(sd$dom))
  expect_error(povmap:::.megb_fold_boot_setup(ff, Xs[so, ], rev(so), as.character(sd$dom[so]),
                                              NULL, NULL, NULL, NULL, use_leaf = FALSE), "out-of-fold")
})

test_that("megb's default bootstrap refreshes every booster with the survey weights and re-estimates by ML from fold models", {
  skip_if_not_installed("xgboost")
  M <- .bf_data(); rec <- new.env(); rec$w_null <- 0L; rec$calls <- 0L; rec$reml <- logical(0)
  orig_r <- povmap:::.megb_refresh; orig_l <- povmap:::.megb_oof_re_fit
  local_mocked_bindings(
    .megb_refresh = function(bst, X, y, w, nr) { rec$calls <- rec$calls + 1L; if (is.null(w)) rec$w_null <- rec$w_null + 1L; orig_r(bst, X, y, w, nr) },
    .megb_oof_re_fit = function(...) { l <- orig_l(...); rec$reml <- c(rec$reml, lme4::isREML(l)); l },
    .package = "povmap")
  f <- .bf_fit(M, B = 4)
  expect_equal(rec$calls, 4L * (1L + 5L))           # per replicate: the full booster and five fold models
  expect_equal(rec$w_null, 0L)
  expect_equal(rec$reml, rep(FALSE, 4L))
  expect_true(all(is.finite(f$var$Mean)))
  expect_true(all(f$CI$Lower <= f$CI$Upper))
})

test_that("the default bootstrap leaves the point estimates unchanged and also runs with row folds", {
  skip_if_not_installed("xgboost")
  M  <- .bf_data()
  f0 <- .bf_fit(M, mse = FALSE); f1 <- .bf_fit(M, B = 3)
  expect_identical(f1$ind$Mean, f0$ind$Mean)
  fr <- .bf_fit(M, B = 3, cv_folds = "rows")
  expect_true(all(is.finite(fr$var$Mean)))
})

test_that("cross-fitted prediction now has a bootstrap: intervals around the cross-fitted estimates", {
  skip_if_not_installed("xgboost")
  M <- .bf_data()
  f <- .bf_fit(M, B = 4, predict_sampled = "crossfit")
  expect_identical(f$CI$Domain, f$ind$Domain)
  expect_true(all(f$CI$Lower <= f$CI$Upper))
  expect_true(all(is.finite(f$var$Mean)))
  expect_error(.bf_fit(M, B = 2, predict_sampled = "crossfit", bootstrap_refit = "lmm_only"), "lmm_only")
})
