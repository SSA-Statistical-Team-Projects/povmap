# megb() defaults from povmap 1.0.1: ten domain folds, cross-fitted prediction for
# sampled domains, the fold-model average for unsampled domains, and bootstraps
# that predict their replicates the same way. predict_unsampled = "foldmean" is
# also an option in its own right.
.df_data <- function(D = 30, n_pop = 20, n_smp = 8, seed = 7) {
  set.seed(seed)
  pop <- data.frame(dom = factor(rep(sprintf("d%02d", 1:D), each = n_pop)),
                    x1 = rnorm(D * n_pop), x2 = runif(D * n_pop), w = rpois(D * n_pop, 20) + 1)
  pop$xm <- rnorm(D)[as.integer(pop$dom)]
  pop$y <- 0.5 * pop$x1 + 0.4 * pop$xm + rnorm(D, 0, 0.3)[as.integer(pop$dom)] + rnorm(nrow(pop))
  smp <- do.call(rbind, lapply(split(pop, pop$dom)[1:25], function(d) d[sample(nrow(d), n_smp), ]))
  smp$wt <- runif(nrow(smp), 0.3, 3)
  list(pop = pop, smp = smp)
}
.df_params <- list(eta = 0.3, max_depth = 3, subsample = 1, nrounds = 100)
.df_fit <- function(M, mse = FALSE, B = 3, ...) suppressWarnings(povmap::megb(
  fixed = y ~ x1 + x2 + xm, smp_data = M$smp, smp_weights = "wt",
  pop_data = M$pop[, c("dom", "x1", "x2", "xm", "w")], pop_weights = "w", domains = "dom",
  transformation = "no", mse = mse, B = B, na.rm = FALSE, seed = 5,
  gradient_params = .df_params, ...))
.Xcols <- c("x1", "x2", "xm")

test_that("the defaults are ten domain folds, cross-fitted prediction and the fold-model average", {
  skip_if_not_installed("xgboost")
  M <- .df_data()
  f <- suppressMessages(.df_fit(M))
  expect_identical(f$predict_sampled, "crossfit")
  expect_identical(f$predict_unsampled, "foldmean")
  expect_identical(f$cv_nfold, 10L)
  expect_length(f$megb_model$fold_models, 10L)
  ## the default estimate is the cross-fitted + fold-model-average one of the same fit
  expect_identical(f$crossfit$ind$Mean, f$crossfit$ind$Mean_foldmean)
})

test_that("where fold-model prediction does not apply, NULL falls back to 'full' with a message; an explicit request is an error", {
  skip_if_not_installed("xgboost")
  M <- .df_data()
  expect_message(f <- .df_fit(M, cv_folds = "rows"), "predict_sampled = \"full\"")
  expect_identical(f$predict_sampled, "full"); expect_identical(f$predict_unsampled, "full")
  expect_message(.df_fit(M, mse = TRUE, B = 2, bootstrap_refit = "lmm_only"), "lmm_only")
  expect_error(suppressMessages(.df_fit(M, cv_folds = "rows", predict_sampled = "crossfit")), "not available")
  expect_error(suppressMessages(.df_fit(M, cv_folds = "rows", predict_unsampled = "foldmean")), "not available")
  ## the settings before 1.0.1 remain available
  f5 <- suppressMessages(.df_fit(M, cv_nfold = 5, predict_sampled = "full", predict_unsampled = "full"))
  expect_null(f5$crossfit); expect_identical(f5$cv_nfold, 5L)
})

test_that("predict_unsampled = 'foldmean': unsampled domains are the mean of the fold models, sampled ones unchanged", {
  skip_if_not_installed("xgboost")
  M <- .df_data()
  f <- suppressMessages(.df_fit(M, predict_sampled = "full", predict_unsampled = "foldmean"))
  cf <- f$crossfit$ind; fm <- f$megb_model$fold_models
  expect_true(any(!cf$sampled))
  ## equal to machine precision (the random effect is added in a different order), not bit for bit
  expect_equal(cf$Mean[cf$sampled], cf$Mean_full[cf$sampled], tolerance = 1e-12)
  for (d in cf$Domain[!cf$sampled]) {
    P <- M$pop[M$pop$dom == d, ]
    g <- rowMeans(sapply(fm, function(b) stats::predict(b, data.matrix(P[, .Xcols]))))
    expect_equal(cf$Mean[cf$Domain == d], sum(g * P$w) / sum(P$w), tolerance = 1e-10)
  }
})

test_that("on unperturbed data the bootstrap's fold-model average for unsampled domains reproduces the point estimate's", {
  skip_if_not_installed("xgboost")
  M <- .df_data()
  smp <- M$smp[, c("dom", .Xcols)]; pop <- M$pop[, c("dom", .Xcols)]
  f <- suppressMessages(povmap:::megb_em(Y = M$smp$y, X = smp[, -1], dom_name = "dom", smp_data = smp,
         pop_data = pop, gradient_params = .df_params, seed = 5, smp_weights_vec = M$smp$wt,
         cv_folds = "domain", predict_sampled = "crossfit", predict_unsampled = "foldmean"))
  mm <- f$megb_model; sd <- f$inp_smp_data$smp_data; Y <- f$inp_smp_data$target_var
  so <- order(as.character(sd$dom)); Xs <- as.matrix(sd[so, f$cov_names_proc, drop = FALSE])
  pp <- f$pop_data_proc; po <- order(as.character(pp$dom)); pp <- pp[po, , drop = FALSE]
  domains <- as.character(unique(pp$dom)); pdi <- match(as.character(pp$dom), domains)
  Xp <- as.matrix(pp[, f$cov_names_proc, drop = FALSE])
  CFB <- povmap:::.megb_fold_boot_setup(povmap:::.megb_fold_fit(mm), Xs, so, as.character(sd$dom[so]), Xp, pdi,
                                        domains, NULL, use_leaf = TRUE, predict = TRUE, predict_unsampled = TRUE)
  st <- povmap:::.megb_fold_boot_step(CFB, Xs, Xp, Y[so], NULL, refresh = FALSE, use_leaf = TRUE,
                                      predict = TRUE, predict_unsampled = TRUE)
  ## megb_em's Indicators are unweighted domain means; unsampled domains carry no random effect
  ## tolerance 1e-7: the leaf aggregation sums leaf values in double precision, xgboost's predict() in single
  ind <- f$Indicators; ind$dom_name <- as.character(ind$dom_name)
  expect_true(length(CFB$posU) > 0)
  expect_equal(st$u_unweighted, ind$Mean[match(domains[CFB$posU], ind$dom_name)], tolerance = 1e-7)
  ## and the sampled domains' booster parts are their (unrefreshed) fold models', as in the point estimate
  u  <- lme4::ranef(mm$effect_model)$dom
  sp <- domains[CFB$pos_all]
  expect_equal(st$unweighted + u[sp, 1], ind$Mean[match(sp, ind$dom_name)], tolerance = 1e-7)
})

test_that("the default bootstrap refreshes all fold models and averages them for unsampled domains", {
  skip_if_not_installed("xgboost")
  M <- .df_data(); rec <- new.env(); rec$pu <- logical(0); rec$refresh <- 0L
  orig_s <- povmap:::.megb_fold_boot_step; orig_r <- povmap:::.megb_refresh
  local_mocked_bindings(
    .megb_fold_boot_step = function(...) { a <- list(...); rec$pu <- c(rec$pu, isTRUE(a$predict_unsampled)); orig_s(...) },
    .megb_refresh = function(...) { rec$refresh <- rec$refresh + 1L; orig_r(...) },
    .package = "povmap")
  f <- suppressMessages(.df_fit(M, mse = TRUE, B = 3))
  expect_equal(rec$pu, rep(TRUE, 3L))
  expect_equal(rec$refresh, 3L * (1L + 10L))           # per replicate: the full booster and ten fold models
  expect_true(all(is.finite(f$var$Mean)))
})

test_that("the full-refit bootstrap predicts each replicate from that replicate's own fold models", {
  skip_if_not_installed("xgboost")
  M <- .df_data(); rec <- new.env(); rec$n <- 0L; rec$args <- list()
  orig <- povmap:::.megb_compose_preds
  local_mocked_bindings(.megb_compose_preds = function(...) {
    rec$n <- rec$n + 1L; a <- list(...); rec$args[[rec$n]] <- c(a$predict_sampled, a$predict_unsampled); orig(...) },
    .package = "povmap")
  f <- suppressMessages(.df_fit(M, mse = TRUE, B = 2, bootstrap_refit = "full", cv_nfold = 5))
  expect_equal(rec$n, 1L + 2L)                          # the point estimate, then each replicate
  expect_true(all(vapply(rec$args, function(z) identical(z, c("crossfit", "foldmean")), logical(1))))
  expect_true(all(is.finite(f$var$Mean)))
})
