# Regression test: the megb bootstrap pairs each area's replicate estimate with the SAME area's replicate truth.
#
# With bootstrap_refit = "full" the replicate estimates used to come out in the level order of the domain factor of the framework's population frame (sampled areas
# first, then the other areas in frame order), while the replicate truths (tau_star) and `domains` are in alphabetical order. .mse_megb_collate() subtracted them row
# by row, so every area was compared with another area and the bootstrap variance, standard errors and intervals were far too large.
#
# Design: one population row per area, 60 areas in scrambled input order, two thirds of them sampled and interleaved with the unsampled ones, a strong area-level covariate
# (area means SD 0.2), near-zero unit noise. For a sampled area the refit reproduces the replicate truth almost exactly (its estimate is the truth plus a small error), so the paired error is tiny; a mispaired area has
# an error of the size of the spread of the area means.
test_that("megb bootstrap replicate estimates are paired by area with the replicate truths", {
  skip_on_cran(); skip_if_not_installed("xgboost"); skip_if_not_installed("lme4")
  skip_if(utils::packageVersion("testthat") < "3.2.0", "needs testthat >= 3.2.0 for local_mocked_bindings")
  set.seed(7)
  D <- 60; i <- 1:D
  x1 <- 0.2 * (i - mean(i)) / sd(i); truth_mean <- 0.4 + 1.0 * x1 + rnorm(D, 0, 0.02)
  pop <- data.frame(dom = sprintf("d%02d", i), x1 = x1, x2 = runif(D), w = 1000 + 50 * i, stringsAsFactors = FALSE)
  sampled <- sprintf("d%02d", i[i %% 3 != 0])
  smp <- do.call(rbind, lapply(sampled, function(d) {
    k <- as.integer(substr(d, 2, 3))
    data.frame(dom = d, x1 = x1[k], x2 = pop$x2[k], y = truth_mean[k] + rnorm(6, 0, 0.01), wt = runif(6, 0.5, 2), stringsAsFactors = FALSE)
  }))
  pop <- pop[sample(D), ]; smp <- smp[sample(nrow(smp)), ]; rownames(pop) <- rownames(smp) <- NULL      # scrambled input order
  gp <- list(eta = 0.3, max_depth = 3, subsample = 1, colsample_bytree = 1, nrounds = 100, nthread = 1)
  old <- Sys.getenv("OMP_NUM_THREADS"); Sys.setenv(OMP_NUM_THREADS = "1"); on.exit(Sys.setenv(OMP_NUM_THREADS = old), add = TRUE)

  orig <- povmap:::mse_megb
  for (refit in c("full", "leaves_only")) {
    cap <- NULL
    testthat::local_mocked_bindings(mse_megb = function(...) { r <- orig(...); cap <<- r; r }, .package = "povmap")
    suppressMessages(suppressWarnings(povmap::megb(
      fixed = y ~ x1 + x2, smp_data = smp, smp_weights = "wt", pop_data = pop, pop_weights = "w", domains = "dom", transformation = "no",
      mse = TRUE, B = 5, bootstrap_refit = refit, bootstrap_cores = 0, na.rm = FALSE, seed = 11, gradient_params = gp, cv_folds = "domain", cv_nfold = 5,
      early_stopping = TRUE, predict_sampled = "crossfit", predict_unsampled = "foldmean", em_maxit = 10, em_tol = 1e-4, mse_type = "var", boot_area_var = "centred")))
    expect_false(is.null(cap), info = refit)
    dm <- as.character(cap$domains)
    expect_identical(dm, sort(dm), info = refit)
    err <- abs(cap$tau_b_orig - cap$tau_star_orig)
    expect_equal(dim(err), c(D, 5), info = refit)
    # sampled areas: the refit recovers the replicate truth, so any area-level mispairing shows as an error of about the SD of the area means (0.2)
    expect_lt(max(err[dm %in% sampled, ]), 0.1, label = paste("max paired error of sampled areas,", refit))
  }
})
