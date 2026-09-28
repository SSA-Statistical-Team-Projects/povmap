# megb's gradient booster must train to the round its cross-validation selected.
# xgboost >= 3 reports that round in cv$early_stop$best_iteration; povmap read
# cv$best_iteration (NULL there) and fell back to the nrounds cap, so every
# booster trained to the cap while the out-of-fold predictions came from fold
# models stopped earlier. This test fails on that behaviour.
test_that("the megb booster stops at the cross-validated best round, not the cap", {
  skip_if_not_installed("xgboost")
  set.seed(3)
  D <- 40; n_pop <- 30; n_smp <- 8
  pop <- data.frame(dom = factor(rep(sprintf("d%02d", 1:D), each = n_pop)),
                    x1 = rnorm(D * n_pop), x2 = runif(D * n_pop), w = rpois(D * n_pop, 20) + 1)
  pop$y <- 0.5 * pop$x1 + rnorm(D, 0, 0.3)[as.integer(pop$dom)] + rnorm(nrow(pop), 0, 1)
  smp <- do.call(rbind, lapply(split(pop, pop$dom)[1:30], function(d) d[sample(nrow(d), n_smp), ]))
  smp$wt <- 1
  cap <- 400L                      # eta 0.3 on noisy data: cross-validation stops long before 400
  f <- suppressMessages(suppressWarnings(povmap::megb(
    fixed = y ~ x1 + x2, smp_data = smp, smp_weights = "wt", pop_data = pop[, c("dom", "x1", "x2", "w")],
    pop_weights = "w", domains = "dom", transformation = "no", mse = FALSE, na.rm = FALSE, seed = 5,
    gradient_params = list(eta = 0.3, max_depth = 4, subsample = 1, nrounds = cap))))
  M <- f$megb_model
  rounds <- xgboost::xgb.get.num.boosted.rounds(M$boosting)
  best   <- which.min(M$eval_log$test_rmse_mean)            # the last EM iteration's cross-validation
  expect_lt(nrow(M$eval_log), cap)                           # early stopping did bind
  expect_lt(rounds, cap)                                     # old behaviour: rounds == cap
  expect_equal(rounds, best)
})
