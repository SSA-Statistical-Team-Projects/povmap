# xgb_cv(configs = ): each fold's model averages the configurations as xgb(configs = ) does, so the
# held-out predictions are those of the averaged estimator; the same set is used in every fold.

make_cv_data <- function(seed = 7) {
  set.seed(seed)
  nd <- 30; ns <- 12
  pop <- expand.grid(sa = seq_len(ns), dom = seq_len(nd))
  pop$sub <- sprintf("%02d_%02d", pop$dom, pop$sa)
  pop$x1 <- rnorm(nrow(pop)); pop$x2 <- runif(nrow(pop))
  u <- rnorm(nd, 0, 0.4)[pop$dom]
  pop$npop <- round(runif(nrow(pop), 50, 500))
  pop$p <- plogis(-1 + 0.8 * pop$x1 + 0.6 * pop$x2 + u)
  pop$dom <- sprintf("%02d", pop$dom)
  sd_ <- sort(sample(unique(pop$dom), 24))
  smp <- do.call(rbind, lapply(sd_, function(d) {
    s <- pop[pop$dom == d, ][sample(ns, 6), ]
    s$y <- rbinom(nrow(s), 25, s$p) / 25
    s$wt <- exp(rnorm(nrow(s), 0, 0.5)) * 100
    s }))
  list(pop = pop[, c("sub", "dom", "x1", "x2", "npop")], smp = smp[, c("sub", "dom", "x1", "x2", "y", "wt")])
}
cv_run <- function(d, ..., cpus = 1) {
  suppressWarnings(suppressMessages(utils::capture.output(res <- povmap::xgb_cv(
    fixed = y ~ x1 + x2, smp_data = d$smp, smp_weights = "wt", pop_data = d$pop, pop_weights = "npop",
    domains = "dom", sub_domains = "sub", transformation = "arcsin", folds = 4, L = 10, seed = 11,
    cpus = cpus, verbose = FALSE, ...))))
  res
}
cfg_a <- data.frame(nrounds = 15, max_depth = 2, eta = 0.3, min_child_weight = 1, lambda = 1,
                    colsample_bytree = 0.6, colsample_bylevel = 1, colsample_bynode = 1, subsample = 0.6)
cfg_b <- transform(cfg_a, nrounds = 40, eta = 0.1, max_depth = 3, min_child_weight = 5)

test_that("xgb_cv with a one-row configs is identical to the scalar call", {
  skip_if_not_installed("xgboost")
  d <- make_cv_data()
  s <- cv_run(d, nrounds = 15, max_depth = 2, eta = 0.3, min_child_weight = 1, lambda = 1,
              colsample_bylevel = 1, colsample_bynode = 1)
  c1 <- cv_run(d, configs = cfg_a)
  expect_identical(c1$domain_results, s$domain_results)
  expect_identical(c1$r2_cv, s$r2_cv)
  expect_equal(c1$configs$weight, 1)
})

test_that("held-out predictions are the weighted average of the configurations' predictions", {
  ## same folds (seed) in every run; without a benchmark xgb()'s Mean is the weighted average of the
  ## configurations' Means, so the averaged CV prediction is the weighted average of the single ones
  skip_if_not_installed("xgboost")
  d <- make_cv_data()
  ca <- cv_run(d, configs = cfg_a); cb <- cv_run(d, configs = cfg_b)
  w <- c(0.3, 0.7)
  cab <- cv_run(d, configs = cbind(rbind(cfg_a, cfg_b), weight = w))
  expect_identical(cab$domain_results$Domain, ca$domain_results$Domain)
  expect_equal(cab$domain_results$Predicted,
               w[1] * ca$domain_results$Predicted + w[2] * cb$domain_results$Predicted, tolerance = 1e-12)
  expect_identical(cab$domain_results$Direct, ca$domain_results$Direct)
  expect_equal(nrow(cab$configs), 2)
})

test_that("xgb_cv with configs gives the same result with 1 and 2 workers", {
  skip_if_not_installed("xgboost")
  skip_on_cran()
  d <- make_cv_data()
  cfgs <- rbind(cfg_a, cfg_b)
  r1 <- cv_run(d, configs = cfgs)
  r2 <- cv_run(d, configs = cfgs, cpus = 2)
  expect_equal(r2$domain_results[order(r2$domain_results$Domain), ],
               r1$domain_results[order(r1$domain_results$Domain), ], ignore_attr = TRUE)
})

test_that("xgb_cv rejects scalar hyperparameters together with configs", {
  skip_if_not_installed("xgboost")
  d <- make_cv_data()
  expect_error(cv_run(d, configs = cfg_a, nrounds = 20), "leave the hyperparameter arguments")
})
