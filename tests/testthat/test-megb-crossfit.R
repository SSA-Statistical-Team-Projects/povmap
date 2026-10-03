# megb(predict_sampled = "crossfit"): each sampled domain's population rows are
# predicted by the fold model of the final EM iteration that held the domain out,
# plus its random effect; unsampled domains keep the full booster. The fit itself
# is unchanged, so predict_sampled = "full" (the default) must reproduce the
# historical estimates and the crossfit fit's Mean_full must equal them.
.cf_data <- function() {
  set.seed(7)
  D <- 30; n_pop <- 20; n_smp <- 8
  pop <- data.frame(dom = factor(rep(sprintf("d%02d", 1:D), each = n_pop)),
                    x1 = rnorm(D * n_pop), x2 = runif(D * n_pop), w = rpois(D * n_pop, 20) + 1)
  pop$xm <- rnorm(D)[as.integer(pop$dom)]          # constant within domain, unique to it
  pop$y <- 0.5 * pop$x1 + 0.4 * pop$xm + rnorm(D, 0, 0.3)[as.integer(pop$dom)] + rnorm(nrow(pop))
  smp <- do.call(rbind, lapply(split(pop, pop$dom)[1:25], function(d) d[sample(nrow(d), n_smp), ]))
  smp$wt <- 1
  list(pop = pop, smp = smp)
}
.cf_fit <- function(M, mse = FALSE, cv_folds = "domain", ...) suppressMessages(suppressWarnings(povmap::megb(
  fixed = y ~ x1 + x2 + xm, smp_data = M$smp, smp_weights = "wt",
  pop_data = M$pop[, c("dom", "x1", "x2", "xm", "w")], pop_weights = "w", domains = "dom",
  transformation = "no", mse = mse, na.rm = FALSE, seed = 5, cv_folds = cv_folds,
  gradient_params = list(eta = 0.3, max_depth = 3, subsample = 1, nrounds = 100), ...)))
.cf_record_cv <- function(env) {
  orig <- xgboost::xgb.cv
  function(...) { env$args <- list(...); orig(...) }   # keeps the last call: the final EM iteration
}

test_that("crossfit: each sampled domain is predicted by the fold model that excluded it", {
  skip_if_not_installed("xgboost")
  M <- .cf_data(); rec <- new.env()
  local_mocked_bindings(xgb.cv = .cf_record_cv(rec), .package = "xgboost")
  orig_tg <- povmap:::train_gbmodel
  local_mocked_bindings(train_gbmodel = function(...) { rec$groups <- as.character(list(...)$groups); orig_tg(...) },
                        .package = "povmap")
  f  <- .cf_fit(M, predict_sampled = "crossfit")
  fm <- f$megb_model$fold_models; fd <- f$megb_model$fold_domains
  expect_length(fm, 5L); expect_length(fd, 5L)
  sampled <- as.character(unique(M$smp$dom))
  ## the folds partition the sampled domains, and fold k's held-out rows (xgb.cv's
  ## test indices for model k, final EM iteration) are exactly its domains' rows
  expect_setequal(unlist(fd), sampled); expect_false(anyDuplicated(unlist(fd)) > 0)
  for (k in 1:5) expect_setequal(which(rec$groups %in% fd[[k]]), rec$args$folds[[k]])
  ## the fold recorded for each sampled domain is the one that held it out
  cf <- f$crossfit$ind
  expect_equal(cf$sampled, cf$Domain %in% sampled)
  for (d in sampled) expect_true(d %in% fd[[cf$fold[cf$Domain == d]]])
  ## sampled rows predicted by their fold model reproduce xgb.cv's out-of-fold predictions
  expect_lt(f$crossfit$oof_max_abs_diff, 1e-6)
  ## each sampled domain's estimate = pop-weighted mean of its fold model's predictions + its random effect
  u  <- lme4::ranef(f$megb_model$effect_model)[[1]]
  for (d in sampled) {
    k <- which(vapply(fd, function(z) d %in% z, logical(1)))
    P <- M$pop[M$pop$dom == d, ]
    g <- stats::predict(fm[[k]], data.matrix(P[, c("x1", "x2", "xm")]))
    expect_equal(cf$Mean[cf$Domain == d], sum((g + u[d, 1]) * P$w) / sum(P$w), tolerance = 1e-10)
  }
  ## and it differs from the full booster's estimate
  expect_true(all(cf$Mean[cf$sampled] != cf$Mean_full[cf$sampled]))
})

test_that("crossfit: unsampled domains keep the full booster; Mean_foldmean averages the five fold models", {
  skip_if_not_installed("xgboost")
  M  <- .cf_data()
  f  <- .cf_fit(M, predict_sampled = "crossfit")
  cf <- f$crossfit$ind; fm <- f$megb_model$fold_models
  expect_true(any(!cf$sampled))
  expect_identical(cf$Mean[!cf$sampled], cf$Mean_full[!cf$sampled])
  expect_identical(cf$Mean_foldmean[cf$sampled], cf$Mean[cf$sampled])
  for (d in cf$Domain[!cf$sampled]) {
    P <- M$pop[M$pop$dom == d, ]
    g <- rowMeans(sapply(fm, function(b) stats::predict(b, data.matrix(P[, c("x1", "x2", "xm")]))))
    expect_equal(cf$Mean_foldmean[cf$Domain == d], sum(g * P$w) / sum(P$w), tolerance = 1e-10)
  }
})

test_that("predict_sampled = 'full' is the default and the fit is the same under both", {
  skip_if_not_installed("xgboost")
  M  <- .cf_data()
  f0 <- .cf_fit(M)
  f1 <- .cf_fit(M, predict_sampled = "crossfit")
  expect_identical(f0$predict_sampled, "full"); expect_null(f0$crossfit)
  expect_null(f0$megb_model$fold_models)
  expect_identical(f1$crossfit$ind$Mean_full, f0$ind$Mean[match(f1$crossfit$ind$Domain, f0$ind$Domain)])
  expect_identical(f1$megb_model$ran_eff_sd, f0$megb_model$ran_eff_sd)
  expect_identical(f1$megb_model$error_sd, f0$megb_model$error_sd)
  expect_identical(f1$megb_model$best_iter, f0$megb_model$best_iter)
  expect_identical(f1$megb_model$iterations_used, f0$megb_model$iterations_used)
})

test_that("crossfit refuses the full-refit bootstrap and row folds", {
  skip_if_not_installed("xgboost")
  M <- .cf_data()
  expect_error(.cf_fit(M, predict_sampled = "crossfit", mse = TRUE, B = 2, bootstrap_refit = "full"), "leaves_only")
  expect_error(.cf_fit(M, predict_sampled = "crossfit", cv_folds = "rows"), "domain")
  expect_error(.cf_fit(M, predict_sampled = "oof"))
})

test_that("cv_nfold = 10: ten folds of whole domains, ten fold models under crossfit", {
  skip_if_not_installed("xgboost")
  M <- .cf_data(); rec <- new.env()
  local_mocked_bindings(xgb.cv = .cf_record_cv(rec), .package = "xgboost")
  f <- .cf_fit(M, predict_sampled = "crossfit", cv_nfold = 10)
  expect_length(rec$args$folds, 10L)
  expect_length(f$megb_model$fold_models, 10L)
  expect_setequal(unlist(f$megb_model$fold_domains), as.character(unique(M$smp$dom)))
  expect_lt(f$crossfit$oof_max_abs_diff, 1e-6)
  expect_error(.cf_fit(M, cv_nfold = 1))
  ## the default is five folds
  .cf_fit(M); expect_length(rec$args$folds, 5L)
})
