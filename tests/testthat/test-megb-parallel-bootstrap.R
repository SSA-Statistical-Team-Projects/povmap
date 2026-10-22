# bootstrap_cores > 1 with bootstrap_refit = "leaves_only" runs the replicates on a PSOCK cluster. The replicate draws no random numbers, so the
# result must be identical (not just close) to the sequential loop for the same seed. Needs the installed package on the workers.
test_that("parallel leaves_only bootstrap is identical to the sequential loop", {
  skip_on_cran(); skip_if_not_installed("xgboost"); skip_if_not_installed("collapse"); skip_if_not_installed("jsonlite")
  skip_if_not("povmap" %in% rownames(utils::installed.packages()))
  set.seed(42)
  D <- 40; n_pop <- 30; n_smp <- 6
  pop <- data.frame(dom = factor(rep(sprintf("d%02d", 1:D), each = n_pop)),
                    x1 = rnorm(D * n_pop), x2 = runif(D * n_pop), w = rpois(D * n_pop, 20) + 1)
  u <- rnorm(D, 0, 0.3)
  pop$y <- 1 + 0.8 * pop$x1 - 0.5 * pop$x2^2 + u[as.integer(pop$dom)] + rnorm(nrow(pop), 0, 0.5)
  smp <- do.call(rbind, lapply(split(pop, pop$dom)[1:30], function(d) d[sample(nrow(d), n_smp), ]))
  smp$wt <- runif(nrow(smp), 0.5, 2)
  gp <- list(eta = 0.1, max_depth = 3, subsample = 0.8, colsample_bytree = 0.9, nrounds = 60, nthread = 1)
  old <- Sys.getenv("OMP_NUM_THREADS"); Sys.setenv(OMP_NUM_THREADS = "1"); on.exit(Sys.setenv(OMP_NUM_THREADS = old), add = TRUE)
  fit <- function(cores, leaf, tr) {
    o <- options(povmap.leaves_only.leaf_aggregate = leaf); on.exit(options(o))
    suppressMessages(suppressWarnings(povmap::megb(
      fixed = y ~ x1 + x2, smp_data = smp, smp_weights = "wt", pop_data = pop[, c("dom", "x1", "x2", "w")],
      pop_weights = "w", domains = "dom", transformation = tr, mse = TRUE, B = 6, bootstrap_refit = "leaves_only",
      bootstrap_cores = cores, na.rm = FALSE, seed = 11, gradient_params = gp)))
  }
  for (cfg in list(list(TRUE, "no"), list(FALSE, "no"))) {   # leaf aggregation and per-cell prediction
    a <- fit(0, cfg[[1]], cfg[[2]]); b <- fit(2, cfg[[1]], cfg[[2]])
    for (nm in c("ind", "MSE", "var", "CI")) expect_identical(a[[nm]], b[[nm]], info = paste(nm, cfg[[1]], cfg[[2]]))
  }
})
