# Configuration averaging: xgb(configs = ), xgb_top_configs() and xgb_tune(grid = ).
#
# Design (sae_runs_20261002/config_avg/PROPOSAL_config_averaging.md): the point estimate is the
# weighted average of the configurations' domain estimates, benchmarked once; in the bootstrap each
# replicate uses one configuration drawn from its own stream, so the replicate streams are unchanged.
# With configs = NULL nothing changes (the real-data reproduction test checks that on the Nigeria
# pov ward model).

make_cfg_data <- function(seed = 7) {
  set.seed(seed)
  nd <- 30; ns <- 12
  pop <- expand.grid(sa = seq_len(ns), dom = seq_len(nd))
  pop$sub <- sprintf("%02d_%02d", pop$dom, pop$sa)
  pop$reg <- sprintf("r%d", (pop$dom - 1) %/% 10 + 1)          # benchmark level: 3 regions
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
  list(pop = pop[, c("sub", "dom", "reg", "x1", "x2", "npop")],
       smp = smp[, c("sub", "dom", "reg", "x1", "x2", "y", "wt")])
}

run_xgb <- function(d, ..., bootstrap = TRUE, B = 6, cpus = 1) {
  suppressWarnings(suppressMessages(utils::capture.output(res <- povmap::xgb(
    fixed = y ~ x1 + x2, smp_data = d$smp, smp_weights = "wt", pop_data = d$pop, pop_weights = "npop",
    domains = "dom", sub_domains = "sub", transformation = "arcsin", bootstrap = bootstrap, B = B,
    L = 10, seed = 11, cpus = cpus, benchmark = "Mean", benchmark_level = "reg", benchmark_type = "ratio",
    perturb_benchmark = TRUE, colsample_bylevel = 1, colsample_bynode = 1, ...))))
  res
}
run_xgb_cfg <- function(d, configs, ..., bootstrap = TRUE, B = 6, cpus = 1) {
  ## with configs the scalar hyperparameters stay at their defaults, so put everything in configs
  suppressWarnings(suppressMessages(utils::capture.output(res <- povmap::xgb(
    fixed = y ~ x1 + x2, smp_data = d$smp, smp_weights = "wt", pop_data = d$pop, pop_weights = "npop",
    domains = "dom", sub_domains = "sub", transformation = "arcsin", bootstrap = bootstrap, B = B,
    L = 10, seed = 11, cpus = cpus, benchmark = "Mean", benchmark_level = "reg", benchmark_type = "ratio",
    perturb_benchmark = TRUE, configs = configs, ...))))
  res
}
cfg_a <- data.frame(nrounds = 15, max_depth = 2, eta = 0.3, min_child_weight = 1, lambda = 1,
                    colsample_bytree = 0.6, colsample_bylevel = 1, colsample_bynode = 1, subsample = 0.6)
cfg_b <- transform(cfg_a, nrounds = 40, eta = 0.1, max_depth = 3, min_child_weight = 5)

test_that("a one-row configs gives exactly the scalar call", {
  skip_if_not_installed("xgboost")
  d <- make_cfg_data()
  s <- run_xgb(d, nrounds = 15, max_depth = 2, eta = 0.3, min_child_weight = 1, lambda = 1)
  c1 <- run_xgb_cfg(d, cfg_a)
  expect_identical(c1$ind, s$ind)
  expect_identical(c1$var, s$var)
  expect_identical(c1$CI, s$CI)
  expect_identical(c1$yhat, s$yhat)
  expect_identical(c1$config_draws, rep(1L, 6))
  expect_equal(c1$configs$weight, 1)
})

test_that("duplicated identical configurations reproduce the single configuration exactly", {
  ## the replicates draw configurations, but all draws give the same inputs and the replicate
  ## streams are unchanged, so ind, var and CI must be identical
  skip_if_not_installed("xgboost")
  d <- make_cfg_data()
  c1 <- run_xgb_cfg(d, cfg_a)
  c3 <- run_xgb_cfg(d, rbind(cfg_a, cfg_a, cfg_a))
  expect_true(length(unique(c3$config_draws)) > 1)
  expect_identical(c3$var, c1$var)                         # the replicates do not involve the point
  ## the point is (x + x + x) / 3, which can differ from x in the last bit, and the intervals are
  ## centred on it
  expect_equal(c3$ind, c1$ind, tolerance = 1e-12)
  expect_equal(c3$CI, c1$CI, tolerance = 1e-12)
})

test_that("the averaged point estimate is the benchmarked weighted mean of the configurations' estimates", {
  skip_if_not_installed("xgboost")
  d <- make_cfg_data()
  ca <- run_xgb_cfg(d, cfg_a, bootstrap = FALSE)
  cb <- run_xgb_cfg(d, cfg_b, bootstrap = FALSE)
  w <- c(0.25, 0.75)
  cab <- run_xgb_cfg(d, cbind(rbind(cfg_a, cfg_b), weight = w), bootstrap = FALSE)
  expect_equal(cab$ind$Mean, w[1] * ca$ind$Mean + w[2] * cb$ind$Mean, tolerance = 1e-12)
  expect_equal(cab$ind$Mean_agg, w[1] * ca$ind$Mean_agg + w[2] * cb$ind$Mean_agg, tolerance = 1e-12)
  expect_equal(cab$ind_by_config$Mean_config1, ca$ind$Mean, tolerance = 1e-12)
  expect_equal(cab$ind_by_config$Mean_config2, cb$ind$Mean, tolerance = 1e-12)
  ## benchmarked once, on the average: state means of Mean_bench hit the targets, and Mean_bench is
  ## add_benchmark() of the averaged Mean
  bm <- povmap:::add_benchmark(x = cab$ind$Mean, benchmark_level = "reg", fwk = cab$framework,
                               fixed = y ~ x1 + x2, benchmark = "Mean", benchmark_type = "ratio")
  expect_equal(cab$ind$Mean_bench, bm, tolerance = 1e-12)
  expect_false(isTRUE(all.equal(cab$ind$Mean, ca$ind$Mean)))
})

test_that("configuration draws are reproducible, follow the weights, and do not depend on cpus", {
  skip_if_not_installed("xgboost")
  skip_on_cran()
  d <- make_cfg_data()
  cfgs <- cbind(rbind(cfg_a, cfg_b), weight = c(0.5, 0.5))
  r1 <- run_xgb_cfg(d, cfgs, B = 8)
  r2 <- run_xgb_cfg(d, cfgs, B = 8)
  expect_identical(r1$config_draws, r2$config_draws)
  expect_identical(r1$var, r2$var)
  r3 <- run_xgb_cfg(d, cfgs, B = 8, cpus = 2)
  expect_identical(r3$config_draws, r1$config_draws)
  expect_identical(r3$var, r1$var)
  expect_identical(r3$CI, r1$CI)
  ## a different config_seed changes the draws only
  r4 <- run_xgb_cfg(d, cfgs, B = 8, config_seed = 99)
  expect_false(identical(r4$config_draws, r1$config_draws))
  expect_identical(r4$ind, r1$ind)
  ## the draws follow the weights
  set.seed(1)
  r5 <- run_xgb_cfg(d, cbind(rbind(cfg_a, cfg_b), weight = c(1, 0)), B = 8)
  expect_true(all(r5$config_draws == 1L))
})

test_that("each replicate uses its drawn configuration's fit and hyperparameters together", {
  ## weights (1, 0): every replicate draws a, so everything equals the single-configuration run of a.
  ## weights (0, 1): every replicate draws b, whose fit generates the population and whose
  ## hyperparameters refit it, so the point and the unbenchmarked bootstrap equal the run of b.
  ## (The benchmark-target perturbations follow the RNG state after configuration 1's point fit, so
  ## the benchmarked columns of the (0, 1) run are not compared.)
  skip_if_not_installed("xgboost")
  d <- make_cfg_data()
  ra <- run_xgb_cfg(d, cfg_a); rb <- run_xgb_cfg(d, cfg_b)
  r10 <- run_xgb_cfg(d, cbind(rbind(cfg_a, cfg_b), weight = c(1, 0)))
  expect_true(all(r10$config_draws == 1L))
  expect_identical(r10$ind, ra$ind); expect_identical(r10$var, ra$var); expect_identical(r10$CI, ra$CI)
  r01 <- run_xgb_cfg(d, cbind(rbind(cfg_a, cfg_b), weight = c(0, 1)))
  expect_true(all(r01$config_draws == 2L))
  expect_identical(r01$ind$Mean, rb$ind$Mean)
  expect_identical(r01$var$Mean, rb$var$Mean)
  expect_identical(r01$CI[, c("Domain", "Lower", "Upper")], rb$CI[, c("Domain", "Lower", "Upper")])
  expect_false(isTRUE(all.equal(ra$var$Mean, rb$var$Mean)))
})

test_that("configs rejects scalar hyperparameters and malformed input", {
  skip_if_not_installed("xgboost")
  d <- make_cfg_data()
  expect_error(run_xgb_cfg(d, cfg_a, nrounds = 20, bootstrap = FALSE), "leave the hyperparameter arguments")
  expect_error(run_xgb_cfg(d, transform(cfg_a, foo = 1), bootstrap = FALSE), "unknown columns")
  expect_error(run_xgb_cfg(d, cfg_a[0, ], bootstrap = FALSE), "one row per configuration")
  expect_error(run_xgb_cfg(d, cbind(cfg_a, weight = -1), bootstrap = FALSE), "weight")
})

## ---- xgb_top_configs ---------------------------------------------------------------------------
fake_tune <- function(cv, grid = NULL) {
  n <- nrow(cv)
  if (is.null(grid)) grid <- data.frame(nround = seq(10, by = 10, length.out = n), max_depth = 2,
                                        colsample_bytree = 1, colsample_bylevel = 1, colsample_bynode = 1,
                                        subsample = 1, min_child_weight = 1, eta = 0.3, gamma = 0,
                                        max_delta_step = 0, lambda = 1, alpha = 0)
  grid$mse <- rowMeans(cv); grid$rank <- rank(grid$mse, ties.method = "first")
  list(candidates = grid, cv_by_fold = cv)
}

test_that("xgb_top_configs: paired one-SE set, floor, cap, k_best and one_se", {
  ## 6 configurations x 4 folds. Config 1 is best; 2 differs by a constant 0.05 (paired SE 0: out);
  ## 3 differs by noise around 0 (in); 4 is far worse; 5 and 6 slightly worse with noise.
  cv <- rbind(c(1.00, 1.20, 0.90, 1.10),
              c(1.05, 1.25, 0.95, 1.15),
              c(1.02, 1.18, 0.93, 1.09),
              c(2.00, 2.20, 1.90, 2.10),
              c(1.03, 1.30, 0.88, 1.12),
              c(1.10, 1.21, 0.97, 1.16))
  t1 <- fake_tune(cv)
  p <- xgb_top_configs(t1, floor = 1, cap = 8)
  ## paired differences from config 1
  d <- sweep(cv, 2, cv[1, ]); keep <- rowMeans(d) <= apply(d, 1, sd) / 2; keep[1] <- TRUE
  expect_setequal(p$nrounds, t1$candidates$nround[which(keep)])
  expect_false(20 %in% p$nrounds)                         # constant-offset config excluded
  expect_equal(p$nrounds[1], 10)                           # best first
  expect_equal(sum(p$weight), 1)
  expect_true(all(c("nrounds", "max_depth", "eta", "weight", "cv_score", "rank") %in% names(p)))
  ## floor fills with the next best by score; cap keeps the best by score
  p3 <- xgb_top_configs(t1, floor = 4, cap = 8)
  expect_equal(nrow(p3), max(4, sum(keep)))
  expect_equal(p3$cv_score, sort(p3$cv_score))
  p2 <- xgb_top_configs(t1, floor = 1, cap = 2)
  expect_equal(nrow(p2), 2); expect_equal(p2$nrounds[1], 10)
  expect_equal(xgb_top_configs(t1, rule = "k_best", k = 3, floor = 1)$nrounds, t1$candidates$nround[order(rowMeans(cv))[1:3]])
  ## one_se: score within best + SE(best)
  se_b <- sd(cv[1, ]) / 2
  expect_setequal(xgb_top_configs(t1, rule = "one_se", floor = 1)$nrounds, t1$candidates$nround[rowMeans(cv) <= mean(cv[1, ]) + se_b])
})

test_that("xgb_top_configs scores over several fold splits and checks the grids match", {
  set.seed(3)
  cvs <- lapply(1:3, function(s) matrix(rep(c(1, 1.01, 1.2, 1.02), 5) + rnorm(20, 0, 0.02), 4))
  ts <- lapply(cvs, fake_tune)
  p <- xgb_top_configs(ts, floor = 1)
  score <- Reduce(`+`, lapply(cvs, rowMeans)) / 3
  expect_equal(p$cv_score[1], min(score))
  expect_equal(p$cv_score, sort(score)[seq_len(nrow(p))])
  ## the paired SE averages the per-split variances (not divided by the number of splits)
  b <- which.min(score)
  dd <- Reduce(`+`, lapply(cvs, function(cf) rowMeans(sweep(cf, 2, cf[b, ])))) / 3
  se <- sqrt(Reduce(`+`, lapply(cvs, function(cf) apply(sweep(cf, 2, cf[b, ]), 1, var) / 5)) / 3)
  expect_setequal(p$nrounds, ts[[1]]$candidates$nround[dd <= se | seq_along(dd) == b])
  bad <- ts; bad[[2]]$candidates$max_depth <- 3
  expect_error(xgb_top_configs(bad), "same grid")
  expect_error(xgb_top_configs(list(list(candidates = NULL))), "keep_candidates")
})

## ---- xgb_tune(grid = ) -----------------------------------------------------------------------------
test_that("xgb_tune(grid = ) scores exactly the configurations of the equivalent full cross", {
  skip_if_not_installed("xgboost")
  d <- make_cfg_data()$smp
  args <- list(y ~ x1 + x2, smp_data = d, smp_weights = "wt", domains = "dom", folds = 3, seed = 5L,
               cpus = 1, verbose = FALSE, colsample_bytree = 1, subsample = 1)
  full <- suppressWarnings(do.call(xgb_tune, c(args, list(nround = c(5, 10), max_depth = c(1, 2), min_child_weight = 1, lambda = 1))))
  g <- expand.grid(nround = c(5, 10), max_depth = c(1, 2))
  viagrid <- suppressWarnings(do.call(xgb_tune, c(args, list(grid = g, min_child_weight = 1, lambda = 1))))
  expect_identical(viagrid$cv_by_fold, full$cv_by_fold)
  expect_identical(viagrid$candidates$mse, full$candidates$mse)
  ## a grid that is not a full cross: more rounds only at the smaller eta
  g2 <- rbind(data.frame(eta = 0.3, nround = 5), data.frame(eta = 0.1, nround = c(5, 15)))
  part <- suppressWarnings(do.call(xgb_tune, c(args, list(grid = g2, max_depth = 2, min_child_weight = 1, lambda = 1))))
  expect_equal(nrow(part$candidates), 3)
  expect_equal(part$candidates$eta, c(0.3, 0.1, 0.1))
  expect_equal(part$candidates$max_depth, rep(2, 3))       # missing columns take the argument's first value
  expect_error(suppressWarnings(do.call(xgb_tune, c(args, list(grid = data.frame(foo = 1))))), "not recognised")
  ## and the selected set can be passed straight to xgb()
  top <- xgb_top_configs(part, floor = 2)
  expect_true(all(c("nrounds", "eta", "weight") %in% names(top)))
})
