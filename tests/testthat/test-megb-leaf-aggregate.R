# options(povmap.leaves_only.leaf_aggregate = TRUE): the leaves_only bootstrap
# computes its domain means from refreshed leaf values instead of predicting every
# population row. The point estimate must be unchanged; the bootstrap outputs may
# differ only by floating-point accumulation (xgboost predicts in single precision).
test_that("leaf-level aggregation reproduces the leaves_only bootstrap", {
  skip_if_not_installed("xgboost"); skip_if_not_installed("collapse"); skip_if_not_installed("jsonlite")
  set.seed(42)
  D <- 40; n_pop <- 30; n_smp <- 6
  pop <- data.frame(dom = factor(rep(sprintf("d%02d", 1:D), each = n_pop)),
                    x1 = rnorm(D * n_pop), x2 = runif(D * n_pop), w = rpois(D * n_pop, 20) + 1)
  u <- rnorm(D, 0, 0.3)
  pop$y <- 1 + 0.8 * pop$x1 - 0.5 * pop$x2^2 + u[as.integer(pop$dom)] + rnorm(nrow(pop), 0, 0.5)
  smp <- do.call(rbind, lapply(split(pop, pop$dom)[1:30], function(d) d[sample(nrow(d), n_smp), ]))
  smp$wt <- runif(nrow(smp), 0.5, 2)
  gp <- list(eta = 0.1, max_depth = 3, subsample = 0.8, colsample_bytree = 0.9, nrounds = 60)
  fit <- function(leaf, pw = "w") {
    old <- options(povmap.leaves_only.leaf_aggregate = leaf); on.exit(options(old))
    suppressMessages(suppressWarnings(povmap::megb(
      fixed = y ~ x1 + x2, smp_data = smp, smp_weights = "wt", pop_data = pop[, c("dom", "x1", "x2", "w")],
      pop_weights = pw, domains = "dom", transformation = "no", mse = TRUE, B = 8,
      na.rm = FALSE, seed = 11, gradient_params = gp)))
  }
  for (pw in list("w", NULL)) {                 # population-weighted and unweighted aggregation
    a <- fit(FALSE, pw); b <- fit(TRUE, pw)
    expect_identical(a$ind, b$ind)              # the point estimate does not involve the bootstrap
    expect_equal(a$MSE, b$MSE, tolerance = 1e-5)
    expect_equal(a$var, b$var, tolerance = 1e-5)
    expect_equal(a$CI, b$CI, tolerance = 1e-5)
  }
})
