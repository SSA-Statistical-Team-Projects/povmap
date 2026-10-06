# megb's bootstrap world under predict_sampled = "crossfit" (5 Oct 2026):
#  * the bootstrap outcome y* is centred on the cross-fitted fit, the fold models'
#    out-of-fold predictions, not on the full booster's in-sample fit;
#  * the residual pool is taken around that cross-fitted fit plus the random effect;
#  * the point estimates do not change;
#  * the full-refit bootstrap's replicate fits use the point fit's EM settings
#    (25 iterations; they were capped at 10).
.bc_data <- function(D = 30, n_pop = 20, n_smp = 8, seed = 7) {
  set.seed(seed)
  pop <- data.frame(dom = factor(rep(sprintf("d%02d", 1:D), each = n_pop)),
                    x1 = rnorm(D * n_pop), x2 = runif(D * n_pop), w = rpois(D * n_pop, 20) + 1)
  pop$xm <- rnorm(D)[as.integer(pop$dom)]
  pop$y <- 0.5 * pop$x1 + 0.4 * pop$xm + rnorm(D, 0, 0.3)[as.integer(pop$dom)] + rnorm(nrow(pop))
  smp <- do.call(rbind, lapply(split(pop, pop$dom)[1:25], function(d) d[sample(nrow(d), n_smp), ]))
  smp$wt <- runif(nrow(smp), 0.3, 3)
  list(pop = pop, smp = smp)
}
.bc_fit <- function(M, mse = TRUE, B = 3, predict_sampled = "crossfit", predict_unsampled = "foldmean",
                    bootstrap_refit = "leaves_only")
  suppressMessages(suppressWarnings(povmap::megb(
  fixed = y ~ x1 + x2 + xm, smp_data = M$smp, smp_weights = "wt",
  pop_data = M$pop[, c("dom", "x1", "x2", "xm", "w")], pop_weights = "w", domains = "dom",
  transformation = "no", mse = mse, B = B, na.rm = FALSE, seed = 5, cv_nfold = 5,
  gradient_params = list(eta = 0.3, max_depth = 3, subsample = 1, nrounds = 100),
  predict_sampled = predict_sampled, predict_unsampled = predict_unsampled,
  bootstrap_refit = bootstrap_refit)))
## each sample row's out-of-fold prediction: the fold model that held its domain out
.bc_oof <- function(model, X, dom) {
  k <- vapply(as.character(dom), function(d) which(vapply(model$fold_domains, function(z) d %in% z, logical(1))), integer(1))
  p <- numeric(nrow(X)); for (j in unique(k)) p[k == j] <- stats::predict(model$fold_models[[j]], X[k == j, , drop = FALSE]); p
}
## records the arguments of povmap-internal function `fn` (all calls), then calls it, or `replace` if given
.bc_record <- function(env, fn, replace = NULL) {
  orig <- utils::getFromNamespace(fn, "povmap")
  function(...) { env[[fn]] <- c(env[[fn]], list(list(...))); if (is.null(replace)) orig(...) else replace(...) }
}

test_that("crossfit: the bootstrap outcome is centred on the out-of-fold predictions", {
  skip_if_not_installed("xgboost")
  M <- .bc_data(); rec <- new.env()
  ## zero residual and random-effect draws, so that y* is exactly the centre of the bootstrap world
  zero_rc <- function(Y, smp_data, unit_pred_smp, ...) {
    smp_data$gb_res <- 0
    list(gb_res = rep(0, length(Y)), ran_effs = rep(0, length(unique(smp_data[[list(...)$dom_name]]))),
         smp_data = smp_data, area_w = NULL)
  }
  ## y* and the sample covariates are read where each replicate's leaf refresh receives them (.megb_leaf_em(CFB, X_smp, X_pop, y, ...))
  local_mocked_bindings(ran_comp = .bc_record(rec, "ran_comp", zero_rc),
                        .megb_leaf_em = .bc_record(rec, ".megb_leaf_em"), .package = "povmap")
  f <- try(.bc_fit(M, B = 2), silent = TRUE)   # replicates on a zero-noise world may fail; only y* is checked
  r <- rec$.megb_leaf_em[[1]]; X <- r[[2]]; y_star <- as.numeric(r[[4]])
  dom <- rec$ran_comp[[1]]$smp_data[[rec$ran_comp[[1]]$dom_name]]
  model <- rec$ran_comp[[1]]$model
  oof  <- .bc_oof(model, X, dom)
  full <- stats::predict(model$boosting, X)
  expect_equal(y_star, oof, tolerance = 1e-6)
  expect_gt(max(abs(y_star - full)), 1e-3)          # the old centre, the full booster's in-sample fit
  ## with predict_sampled = "full" the centre stays the full booster's
  rec2 <- new.env()
  local_mocked_bindings(ran_comp = .bc_record(rec2, "ran_comp", zero_rc),
                        .megb_leaf_em = .bc_record(rec2, ".megb_leaf_em"), .package = "povmap")
  f2 <- try(.bc_fit(M, B = 2, predict_sampled = "full", predict_unsampled = "full"), silent = TRUE)
  r2 <- rec2$.megb_leaf_em[[1]]
  expect_equal(as.numeric(r2[[4]]), unname(as.numeric(stats::predict(rec2$ran_comp[[1]]$model$boosting, r2[[2]]))), tolerance = 1e-6)
})

test_that("crossfit: the residual pool is taken around the cross-fitted fit plus the random effect", {
  skip_if_not_installed("xgboost")
  M <- .bc_data(); rec <- new.env()
  local_mocked_bindings(ran_comp = .bc_record(rec, "ran_comp"), .megb_leaf_em = .bc_record(rec, ".megb_leaf_em"),
                        .package = "povmap")
  f <- .bc_fit(M, B = 2)
  a <- rec$ran_comp[[1]]; X <- rec$.megb_leaf_em[[1]][[2]]
  dom <- a$smp_data[[a$dom_name]]; em <- a$model$effect_model
  re  <- stats::predict(em, a$smp_data, allow.new.levels = TRUE) - lme4::fixef(em)
  oof <- .bc_oof(a$model, X, dom)
  ## the pool is Y minus what ran_comp is given: the out-of-fold fit plus the point fit's random effect
  expect_equal(as.numeric(a$Y - a$unit_pred_smp), as.numeric(a$Y - oof - re), tolerance = 1e-6)
  expect_gt(max(abs(a$unit_pred_smp - (stats::predict(a$model$boosting, X) + re))), 1e-3)
})

test_that("the point estimates do not depend on the bootstrap or its centring", {
  skip_if_not_installed("xgboost")
  M <- .bc_data()
  p0 <- .bc_fit(M, mse = FALSE)$ind
  for (br in c("leaves_only", "full")) {
    f <- .bc_fit(M, B = 2, bootstrap_refit = br)
    expect_identical(f$ind, p0)
  }
  expect_identical(.bc_fit(M, B = 2, predict_sampled = "full", predict_unsampled = "full")$ind,
                   .bc_fit(M, mse = FALSE, predict_sampled = "full", predict_unsampled = "full")$ind)
})

test_that("the full-refit replicates use the point fit's EM settings", {
  skip_if_not_installed("xgboost")
  M <- .bc_data(); rec <- new.env()
  local_mocked_bindings(em_gb_lmm = .bc_record(rec, "em_gb_lmm"), .package = "povmap")
  f <- .bc_fit(M, B = 2, bootstrap_refit = "full")
  calls <- rec$em_gb_lmm
  expect_length(calls, 3L)                          # the point fit, then one per replicate
  it  <- vapply(calls, function(z) as.numeric(z$max_iterations), numeric(1))
  tol <- vapply(calls, function(z) as.numeric(z$error_tolerance), numeric(1))
  expect_equal(it, rep(25, 3)); expect_equal(tol, rep(1e-4, 3))
  expect_identical(povmap:::.megb_em_max_iterations, 25L)
})
