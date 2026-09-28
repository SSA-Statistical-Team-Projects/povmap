# megb trains the full-data xgboost booster only after the EM converges
# (options(megb.defer_full_fit), default TRUE). The estimates, the fitted model
# and the caller's RNG stream must be identical to the per-iteration fit.
test_that("deferring the full-data booster to the last EM iteration changes nothing", {
  skip_if_not_installed("xgboost")
  set.seed(42)
  D <- 40; n_pop <- 30; n_smp <- 6
  pop <- data.frame(dom = factor(rep(sprintf("d%02d", 1:D), each = n_pop)),
                    x1 = rnorm(D * n_pop), x2 = runif(D * n_pop), w = rpois(D * n_pop, 20) + 1)
  u <- rnorm(D, 0, 0.3)
  pop$y <- 1 + 0.8 * pop$x1 - 0.5 * pop$x2^2 + u[as.integer(pop$dom)] + rnorm(nrow(pop), 0, 0.5)
  smp <- do.call(rbind, lapply(split(pop, pop$dom)[1:30], function(d) d[sample(nrow(d), n_smp), ]))
  smp$wt <- runif(nrow(smp), 0.5, 2)
  gp <- list(eta = 0.1, max_depth = 3, subsample = 0.8, colsample_bytree = 0.9, nrounds = 60)
  fit <- function(defer) {
    old <- options(megb.defer_full_fit = defer); on.exit(options(old))
    set.seed(7)
    f <- suppressMessages(suppressWarnings(povmap::megb(
      fixed = y ~ x1 + x2, smp_data = smp, smp_weights = "wt", pop_data = pop[, c("dom", "x1", "x2", "w")],
      pop_weights = "w", domains = "dom", transformation = "no", mse = TRUE, B = 5,
      na.rm = FALSE, seed = 11, gradient_params = gp)))
    list(f = f, after = runif(3))
  }
  a <- fit(FALSE); b <- fit(TRUE)
  expect_identical(a$f$ind, b$f$ind)
  expect_identical(a$f$MSE, b$f$MSE)
  expect_identical(a$f$megb_model$iterations_used, b$f$megb_model$iterations_used)
  expect_identical(a$f$megb_model$error_sd, b$f$megb_model$error_sd)
  expect_identical(a$f$megb_model$ran_eff_sd, b$f$megb_model$ran_eff_sd)
  expect_identical(a$after, b$after)            # the caller's RNG stream is unchanged
})
