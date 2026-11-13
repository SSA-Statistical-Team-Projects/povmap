# Regression test: the megb bootstrap pairs each area's replicate estimate with the SAME area's replicate truth.
#
# With bootstrap_refit = "full" the replicate estimates used to come out in the level order of the domain factor of the framework's population frame (sampled areas
# first, then the other areas in frame order), while the replicate truths (tau_star) and `domains` are in alphabetical order. .mse_megb_collate() subtracted them row
# by row, so every area was compared with another area and the bootstrap variance, standard errors and intervals were far too large.
#
# Design: one population row per area, 60 areas in scrambled input order, two thirds of them sampled and interleaved with the unsampled ones, a strong area-level covariate
# (area means SD 0.2), near-zero unit noise. For a sampled area the refit reproduces the replicate truth almost exactly (its estimate is the truth plus a small error), so the
# paired error is small; a mispaired area has an error of the size of the spread of the area means (0.2).
# Two scenarios per refit type: (a) no transformation, no benchmark; (b) arcsine transformation (per-cell back-transform, the mean_boot_orig branch) and ratio benchmarking
# to four groups with each replicate's own truth as the target (bench_target = "truth"; the benchmark_fn branch), checking tau_b_orig and tau_b_bench. With a fixed target the benchmark factor itself moves the sampled areas by about 0.17, which would hide a pairing error.
test_that("megb bootstrap replicate estimates are paired by area with the replicate truths", {
  skip_on_cran(); skip_if_not_installed("xgboost"); skip_if_not_installed("lme4")
  skip_if(utils::packageVersion("testthat") < "3.2.0", "needs testthat >= 3.2.0 for local_mocked_bindings")
  set.seed(7)
  D <- 60; i <- 1:D
  x1 <- 0.2 * (i - mean(i)) / sd(i); truth_mean <- 0.4 + 1.0 * x1 + rnorm(D, 0, 0.02)
  pop <- data.frame(dom = sprintf("d%02d", i), x1 = x1, x2 = runif(D), w = 1000 + 50 * i, grp = paste0("g", (i - 1) %/% 15 + 1), stringsAsFactors = FALSE)
  sampled <- sprintf("d%02d", i[i %% 3 != 0])
  smp <- do.call(rbind, lapply(sampled, function(d) {
    k <- as.integer(substr(d, 2, 3))
    data.frame(dom = d, x1 = x1[k], x2 = pop$x2[k], grp = pop$grp[k], y = truth_mean[k] + rnorm(6, 0, 0.01), wt = runif(6, 0.5, 2), stringsAsFactors = FALSE)
  }))
  bm <- vapply(split(seq_len(D), pop$grp), function(ix) stats::weighted.mean(truth_mean[ix], pop$w[ix]), 0)   # the true group means as the benchmark target
  pop <- pop[sample(D), ]; smp <- smp[sample(nrow(smp)), ]; rownames(pop) <- rownames(smp) <- NULL      # scrambled input order
  gp <- list(eta = 0.3, max_depth = 3, subsample = 1, colsample_bytree = 1, nrounds = 100, nthread = 1)
  old <- Sys.getenv("OMP_NUM_THREADS"); Sys.setenv(OMP_NUM_THREADS = "1"); on.exit(Sys.setenv(OMP_NUM_THREADS = old), add = TRUE)

  orig <- povmap:::mse_megb
  for (refit in c("full", "leaves_only", "lmm_only")) for (scen in c("a", "b")) {
    cap <- NULL
    testthat::local_mocked_bindings(mse_megb = function(...) { r <- orig(...); cap <<- r; r }, .package = "povmap")
    pred <- if (refit == "lmm_only") "full" else NULL            # lmm_only has no fold-model bootstrap, so it needs the non-fold predictions
    args <- list(fixed = y ~ x1 + x2, smp_data = smp, smp_weights = "wt", pop_data = pop, pop_weights = "w", domains = "dom", mse = TRUE, B = 5, bootstrap_refit = refit,
                 bootstrap_cores = 0, na.rm = FALSE, seed = 11, gradient_params = gp, cv_folds = "domain", cv_nfold = 5, early_stopping = TRUE, em_maxit = 10, em_tol = 1e-4,
                 mse_type = "var", boot_area_var = "centred")
    if (!is.null(pred)) args <- c(args, list(predict_sampled = pred, predict_unsampled = pred)) else args <- c(args, list(predict_sampled = "crossfit", predict_unsampled = "foldmean"))
    args <- c(args, if (scen == "a") list(transformation = "no") else list(transformation = "arcsin", arcsin_backtransform = "plain", benchmark = bm, benchmark_type = "ratio", benchmark_level = "grp", bench_target = "truth", perturb_benchmark = FALSE, benchmark_target_se = stats::setNames(rep(0.001, length(bm)), names(bm))))
    suppressMessages(suppressWarnings(do.call(povmap::megb, args)))
    lab <- paste(refit, scen)
    expect_false(is.null(cap), info = lab)
    dm <- as.character(cap$domains); expect_identical(dm, sort(dm), info = lab)
    s <- dm %in% sampled
    err <- abs(cap$tau_b_orig - cap$tau_star_orig)
    expect_equal(dim(err), c(D, 5), info = lab)
    expect_lt(mean(err[s, ]), 0.05, label = paste("mean paired error of sampled areas, tau_b_orig,", lab)); expect_lt(max(err[s, ]), 0.15, label = paste("max paired error of sampled areas, tau_b_orig,", lab))
    if (scen == "b") {
      expect_false(is.null(cap$tau_b_bench), info = lab)
      eb <- abs(cap$tau_b_bench - cap$tau_star_orig[, cap$bench_keep, drop = FALSE])
      expect_lt(mean(eb[s, ]), 0.05, label = paste("mean paired error of sampled areas, tau_b_bench,", lab)); expect_lt(max(eb[s, ]), 0.15, label = paste("max paired error of sampled areas, tau_b_bench,", lab))
    }
  }
})
