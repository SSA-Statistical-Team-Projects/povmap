# Out-of-fold residual pools in xgb(): area_pool, sub_pool, area_deflate, sub_deflate, pool_tune.
#
# Everything is opt-in. The pools replace the in-sample residual pools only for the bootstrap population, after the point
# estimate, with the random-number state saved and restored: point estimates do not depend on the pool, and an existing call
# is unchanged. (Bit-identity of the full estimates object against the pre-change build is checked outside the suite by
# the comparison recorded in NEWS; here the default call is checked against an explicit default call.)

make_pool_data <- function(seed = 7) {
  set.seed(seed)
  nd <- 30; ns <- 12
  pop <- expand.grid(sa = seq_len(ns), dom = seq_len(nd))
  pop$sub <- sprintf("%02d_%02d", pop$dom, pop$sa)
  pop$reg <- sprintf("r%d", (pop$dom - 1) %/% 10 + 1)
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
    s$vy <- runif(nrow(s), 0.001, 0.01)
    s }))
  rownames(smp) <- NULL
  list(pop = pop[, c("sub", "dom", "reg", "x1", "x2", "npop")],
       smp = smp[, c("sub", "dom", "reg", "x1", "x2", "y", "wt", "vy")])
}
pool_fit <- function(d, ..., bootstrap = FALSE, B = 4, cfg = NULL, quiet = TRUE, hp = list(nrounds = 15, max_depth = 2, eta = 0.3, min_child_weight = 1,
                     lambda = 1, colsample_bytree = 0.6, colsample_bylevel = 1, colsample_bynode = 1, subsample = 0.6)) {
  a <- list(fixed = y ~ x1 + x2, smp_data = d$smp, smp_weights = "wt", pop_data = d$pop, pop_weights = "npop",
            domains = "dom", sub_domains = "sub", transformation = "arcsin", bootstrap = bootstrap, B = B, L = 10, seed = 11)
  a <- utils::modifyList(c(a, if (is.null(cfg)) hp else list(configs = cfg)), list(...))
  if (quiet) suppressWarnings(suppressMessages(utils::capture.output(res <- do.call(povmap::xgb, a))))
  else suppressMessages(utils::capture.output(res <- do.call(povmap::xgb, a)))
  res
}
cfg_a <- data.frame(nrounds = 15, max_depth = 2, eta = 0.3, min_child_weight = 1, lambda = 1,
                    colsample_bytree = 0.6, colsample_bylevel = 1, colsample_bynode = 1, subsample = 0.6)
cfg_b <- transform(cfg_a, nrounds = 30, eta = 0.1, max_depth = 3, min_child_weight = 3)

test_that("defaults: explicit default arguments reproduce the plain call, and no pools element is added", {
  skip_if_not_installed("xgboost")
  d <- make_pool_data()
  e0 <- pool_fit(d, bootstrap = TRUE)
  e1 <- pool_fit(d, bootstrap = TRUE, area_pool = "in_sample", sub_pool = "in_sample", area_deflate = FALSE, sub_deflate = FALSE)
  ## the recorded call differs by the explicit arguments, and a booster is a pointer (never identical, even for two identical calls)
  e0$out_call <- e1$out_call <- NULL; e0$model <- e1$model <- NULL
  expect_identical(e1, e0)
  expect_null(e0$pools)
})

test_that("fold assignments: area folds by domain, sub-area folds distinct within a domain, RNG state untouched", {
  dom <- rep(sprintf("d%02d", 1:25), each = 4); cell <- paste(dom, rep(1:4, 25), sep = "_")
  set.seed(5); s0 <- .Random.seed
  fa <- povmap:::.xgb_pool_fold_domain(dom, 10, 123)
  fs <- povmap:::.xgb_pool_fold_sub(cell, dom, 10, 123)
  expect_identical(.Random.seed, s0)
  expect_true(all(tapply(fa, dom, function(z) length(unique(z))) == 1L))      # a domain is in one fold
  expect_setequal(unique(fa), 1:10)
  expect_true(all(tapply(fs, dom, function(z) length(unique(z))) == 4L))      # 4 sub-areas, 4 different folds
  ## the assignment is seeded and reproducible; a data.frame assignment is used as given
  expect_identical(fa, povmap:::.xgb_pool_fold_domain(dom, 10, 123))
  as_df <- data.frame(dom = unique(dom), fold = fa[match(unique(dom), dom)])
  expect_identical(povmap:::.xgb_pool_fold_domain(dom, 10, 1, as_df), fa)
  expect_error(povmap:::.xgb_pool_fold_domain(dom, 30, 1), "between 2")
})

test_that("point estimates do not depend on the pool (max abs diff < 1e-10); the intervals do", {
  skip_if_not_installed("xgboost")
  d <- make_pool_data()
  ref <- pool_fit(d, bootstrap = TRUE, benchmark = "Mean", benchmark_level = "reg", benchmark_type = "ratio")
  for (args in list(list(area_pool = "out_of_fold"), list(sub_pool = "out_of_fold"), list(sub_pool = "out_of_fold_domain"),
                    list(area_pool = "out_of_fold", sub_pool = "out_of_fold", area_deflate = TRUE, sub_deflate = TRUE, sub_sampling_var = "vy"),
                    list(area_deflate = TRUE))) {
    e <- do.call(pool_fit, c(list(d, bootstrap = TRUE, benchmark = "Mean", benchmark_level = "reg", benchmark_type = "ratio"), args))
    expect_lt(max(abs(e$ind$Mean - ref$ind$Mean)), 1e-10)
    expect_lt(max(abs(e$ind$Mean_bench - ref$ind$Mean_bench)), 1e-10)
    expect_false(isTRUE(all.equal(e$CI$Lower, ref$CI$Lower)))
  }
})

test_that("out-of-fold area pool equals a hand computation with plain xgboost", {
  skip_if_not_installed("xgboost")
  d <- make_pool_data()
  e <- pool_fit(d, area_pool = "out_of_fold", pool_folds = 5, pool_seed = 3)
  pl <- e$pools$configs[[1]]
  ## folds: the recorded assignment, balanced over the 24 sampled domains
  af <- e$pools$area_folds
  expect_equal(nrow(af), 24L); expect_equal(sort(unique(af$fold)), 1:5)
  smp <- d$smp; fold <- af$fold[match(smp$dom, af$domain)]
  yt <- asin(sqrt(smp$y)); X <- as.matrix(smp[, c("x1", "x2")])
  cellx <- as.matrix(d$pop[match(smp$sub, d$pop$sub), c("x1", "x2")])
  p <- list(max_depth = 2, colsample_bytree = 0.6, colsample_bylevel = 1, colsample_bynode = 1, subsample = 0.6, min_child_weight = 1,
            eta = 0.3, gamma = 0, max_delta_step = 0, lambda = 1, alpha = 0, nthread = 1)
  hat <- numeric(nrow(smp))
  for (k in 1:5) {
    tr <- fold != k; w <- smp$wt[tr]
    w <- w * ave(w, smp$dom[tr], FUN = length) / ave(w, smp$dom[tr], FUN = sum); w <- w / mean(w); w <- w / mean(w)
    set.seed(11)
    fit <- xgboost::xgb.train(xgboost::xgb.DMatrix(X[tr, ], label = yt[tr], weight = w), params = p, nrounds = 15, verbose = 0)
    hat[!tr] <- predict(fit, cellx[!tr, , drop = FALSE])
  }
  r <- yt - hat
  rd <- tapply(r * smp$wt, smp$dom, sum) / tapply(smp$wt, smp$dom, sum)
  expect_equal(unname(pl$resid_domains[match(names(rd), pl$domains)]), as.numeric(rd), tolerance = 1e-10)
  expect_equal(unname(pl$wt_domains[match(names(rd), pl$domains)]), as.numeric(tapply(smp$wt, smp$dom, sum)), tolerance = 1e-12)
  ## the in-sample sub-area pool is untouched by an area-only choice
  e0 <- pool_fit(d, sub_pool = "in_sample", area_deflate = TRUE)
  expect_equal(e$pools$configs[[1]]$resid_sub_domains, e0$pools$configs[[1]]$resid_sub_domains, tolerance = 1e-12)
})

test_that("the sub-area pool: sub-areas of a domain are held out in different folds; residuals equal a hand computation", {
  skip_if_not_installed("xgboost")
  d <- make_pool_data()
  e <- pool_fit(d, sub_pool = "out_of_fold", pool_folds = 6, pool_seed = 4)
  sf <- e$pools$sub_folds
  expect_true(all(tapply(sf$fold, sf$domain, function(z) length(unique(z))) == 6L))      # 6 sub-areas per domain, 6 folds
  smp <- d$smp; fold <- sf$fold[match(smp$sub, sf$sub_domain)]
  yt <- asin(sqrt(smp$y)); X <- as.matrix(smp[, c("x1", "x2")]); cellx <- as.matrix(d$pop[match(smp$sub, d$pop$sub), c("x1", "x2")])
  p <- list(max_depth = 2, colsample_bytree = 0.6, colsample_bylevel = 1, colsample_bynode = 1, subsample = 0.6, min_child_weight = 1,
            eta = 0.3, gamma = 0, max_delta_step = 0, lambda = 1, alpha = 0, nthread = 1)
  hat <- numeric(nrow(smp))
  for (k in 1:6) {
    tr <- fold != k; w <- smp$wt[tr]
    w <- w * ave(w, smp$dom[tr], FUN = length) / ave(w, smp$dom[tr], FUN = sum); w <- w / mean(w); w <- w / mean(w)
    set.seed(11)
    fit <- xgboost::xgb.train(xgboost::xgb.DMatrix(X[tr, ], label = yt[tr], weight = w), params = p, nrounds = 15, verbose = 0)
    hat[!tr] <- predict(fit, cellx[!tr, , drop = FALSE])
  }
  r <- yt - hat
  rs <- r - (ave(r * smp$wt, smp$dom, FUN = sum) / ave(smp$wt, smp$dom, FUN = sum))
  expect_equal(unname(e$pools$configs[[1]]$resid_sub_domains), unname(rs), tolerance = 1e-10)
  ## sub_pool = "out_of_fold_domain" is the sub-area residual of the domain-grouped fits (the DRC scratch-hook construction)
  e2 <- pool_fit(d, sub_pool = "out_of_fold_domain")
  e3 <- pool_fit(d, area_pool = "out_of_fold", sub_pool = "out_of_fold_domain")
  expect_equal(e2$pools$configs[[1]]$resid_sub_domains, e3$pools$configs[[1]]$resid_sub_domains, tolerance = 1e-12)
})

test_that("reproduces the construction of the DRC scratch hook (oof_lib.R: balanced domain folds, set.seed before each fit)", {
  skip_if_not_installed("xgboost")
  d <- make_pool_data()
  e <- pool_fit(d, area_pool = "out_of_fold", sub_pool = "out_of_fold_domain", pool_folds = 10, pool_seed = 123)
  ## oof_pools(): ud <- unique(dom); set.seed(inner_seed); sample(rep(1:inner_folds, length.out = length(ud)))
  smp <- d$smp; ud <- unique(smp$dom); set.seed(123); fdf <- data.frame(Domain = ud, fold = sample(rep(1:10, length.out = length(ud))))
  efold <- fdf$fold[match(smp$dom, fdf$Domain)]
  expect_identical(e$pools$area_folds$fold, fdf$fold)
  yt <- asin(sqrt(smp$y)); X <- as.matrix(smp[, c("x1", "x2")]); Xp <- as.matrix(d$pop[, c("x1", "x2")]); rows <- match(smp$sub, d$pop$sub)
  p <- list(max_depth = 2, colsample_bytree = 0.6, colsample_bylevel = 1, colsample_bynode = 1, subsample = 0.6, min_child_weight = 1,
            eta = 0.3, gamma = 0, max_delta_step = 0, lambda = 1, alpha = 0, nthread = 1)
  hat <- rep(NA_real_, nrow(smp))
  for (k in 1:10) {
    te <- which(efold == k); idx <- which(efold != k); wv <- smp$wt[idx]
    ds <- ave(wv, smp$dom[idx], FUN = sum); dn <- ave(wv, smp$dom[idx], FUN = length); wr <- wv * dn / ds; wr <- wr / mean(wr)
    set.seed(11)
    fit <- xgboost:::xgb.train(data = xgboost:::xgb.DMatrix(data = X[idx, , drop = FALSE], label = yt[idx], weight = wr * 1 / mean(wr * 1)),
                               params = p, nrounds = 15, verbose = 0)
    hat[te] <- predict(fit, Xp[rows[te], , drop = FALSE])
  }
  r <- yt - hat; w <- smp$wt; g <- smp$dom
  rd <- collapse:::fmean(r, g = g, w = w, use.g.names = TRUE); wts <- collapse:::fsum(w, g = g)
  pl <- e$pools$configs[[1]]
  expect_equal(pl$resid_domains, as.numeric(rd), tolerance = 1e-10)
  expect_equal(pl$wt_domains, as.numeric(wts), tolerance = 1e-12)
  expect_equal(pl$resid_sub_domains, r - rd[match(g, names(rd))], tolerance = 1e-10)
})

test_that("deflation: factor, scale conversion, input formats, floor, missing domains", {
  skip_if_not_installed("xgboost")
  d <- make_pool_data()
  base <- pool_fit(d, area_pool = "out_of_fold", area_deflate = FALSE)
  v <- stats::setNames(runif(24, 1e-4, 6e-4), sort(unique(d$smp$dom)))
  e <- pool_fit(d, area_pool = "out_of_fold", area_deflate = TRUE, area_sampling_var = v)
  b <- base$pools$configs[[1]]; z <- e$pools$configs[[1]]
  ## by hand: variance moved to the arcsin scale with the package slope, weighted by the draw weights
  est <- tapply(d$smp$y * d$smp$wt, d$smp$dom, sum) / tapply(d$smp$wt, d$smp$dom, sum)
  vt <- vapply(seq_along(v), function(i) v[[i]] * povmap:::.perturb_slope(est[[names(v)[i]]], sqrt(v[[i]]), "arcsin", povmap:::arcsin_transform)^2, 1)
  names(vt) <- names(v)
  msr <- weighted.mean(b$resid_domains^2, b$wt_domains)
  msv <- weighted.mean(vt[b$domains], b$wt_domains)
  f <- sqrt(max(msr - msv, 0) / msr)
  expect_equal(z$area_deflation$factor, f, tolerance = 1e-10)
  expect_equal(z$resid_domains, b$resid_domains * f, tolerance = 1e-12)
  expect_equal(z$area_deflation$msr, msr, tolerance = 1e-12); expect_equal(z$area_deflation$msv, msv, tolerance = 1e-12)
  ## floored at zero when the sampling variance exceeds the pool
  big <- pool_fit(d, area_pool = "out_of_fold", area_deflate = TRUE, area_sampling_var = v * 1e6)
  expect_equal(big$pools$configs[[1]]$area_deflation$factor, 0); expect_true(all(big$pools$configs[[1]]$resid_domains == 0))
  ## a column of smp_data constant within domain gives the same result as the named vector
  d2 <- d; d2$smp$av <- v[d$smp$dom]
  e2 <- pool_fit(d2, area_pool = "out_of_fold", area_deflate = TRUE, area_sampling_var = "av")
  expect_equal(e2$pools$configs[[1]]$resid_domains, z$resid_domains, tolerance = 1e-12)
  d2$smp$av[1] <- d2$smp$av[1] * 2
  expect_error(pool_fit(d2, area_pool = "out_of_fold", area_deflate = TRUE, area_sampling_var = "av"), "constant within domain")
  ## domains missing from the vector are left out of both means, with a warning
  expect_warning(pool_fit(d, area_pool = "out_of_fold", area_deflate = TRUE, area_sampling_var = v[-(1:3)], quiet = FALSE), "no sampling variance")
  expect_error(pool_fit(d, area_deflate = TRUE, area_sampling_var = c(1, 2)), "named")
  ## default: cluster-robust variance of the domain mean from the sample rows
  hand <- vapply(split(seq_len(nrow(d$smp)), d$smp$dom), function(i) { w <- d$smp$wt[i]; y <- d$smp$y[i]; yb <- sum(w * y) / sum(w)
    length(i) / (length(i) - 1) * sum((w * (y - yb))^2) / sum(w)^2 }, 1)
  expect_equal(povmap:::.xgb_domain_mean_var(d$smp$y, d$smp$wt, d$smp$dom)[names(hand)], hand, tolerance = 1e-12)
  expect_true(is.finite(pool_fit(d, area_deflate = TRUE)$pools$configs[[1]]$area_deflation$factor))
  ## sub pool: variance_y column by default, vector of the same length otherwise
  s1 <- pool_fit(d, sub_deflate = TRUE, variance_y = "vy")$pools$configs[[1]]$sub_deflation
  s2 <- pool_fit(d, sub_deflate = TRUE, sub_sampling_var = d$smp$vy)$pools$configs[[1]]$sub_deflation
  expect_equal(s1$factor, s2$factor, tolerance = 1e-12)
  expect_error(pool_fit(d, sub_deflate = TRUE), "variance_y")
  expect_error(pool_fit(d, sub_deflate = TRUE, sub_sampling_var = 1:3), "one element per row")
  ## in-sample pools can be deflated too
  expect_true(z$area_deflation$factor <= 1)
})

test_that("configuration averaging: one pool per configuration, each as the scalar call with that configuration", {
  skip_if_not_installed("xgboost")
  d <- make_pool_data()
  cfg <- rbind(cfg_a, cfg_b)
  e <- pool_fit(d, cfg = cfg, area_pool = "out_of_fold", sub_pool = "out_of_fold", bootstrap = TRUE)
  expect_length(e$pools$configs, 2L)
  sa <- pool_fit(d, area_pool = "out_of_fold", sub_pool = "out_of_fold")
  sb <- pool_fit(d, area_pool = "out_of_fold", sub_pool = "out_of_fold",
                 hp = list(nrounds = 30, max_depth = 3, eta = 0.1, min_child_weight = 3, lambda = 1, colsample_bytree = 0.6,
                           colsample_bylevel = 1, colsample_bynode = 1, subsample = 0.6))
  expect_equal(e$pools$configs[[1]]$resid_domains, sa$pools$configs[[1]]$resid_domains, tolerance = 1e-12)
  expect_equal(e$pools$configs[[2]]$resid_domains, sb$pools$configs[[1]]$resid_domains, tolerance = 1e-12)
  expect_equal(e$pools$configs[[2]]$resid_sub_domains, sb$pools$configs[[1]]$resid_sub_domains, tolerance = 1e-12)
  expect_false(isTRUE(all.equal(e$pools$configs[[1]]$resid_domains, e$pools$configs[[2]]$resid_domains)))
  ## point estimate of the average unchanged by the pools
  e0 <- pool_fit(d, cfg = cfg, bootstrap = TRUE)
  expect_lt(max(abs(e$ind$Mean - e0$ind$Mean)), 1e-10)
})

test_that("pool_tune: the area pool from the tuning folds equals the pool from freshly fitted identical folds", {
  skip_if_not_installed("xgboost")
  skip_on_cran()
  d <- make_pool_data()
  g1 <- data.frame(nround = 15, max_depth = 2, colsample_bytree = 0.6, colsample_bylevel = 1, colsample_bynode = 1, subsample = 0.6,
                   min_child_weight = 1, eta = 0.3, gamma = 0, max_delta_step = 0, lambda = 1, alpha = 0)
  tune <- suppressWarnings(suppressMessages(povmap::xgb_tune(
    y ~ x1 + x2, smp_data = d$smp, smp_weights = "wt", domains = "dom", transformation = "arcsin", folds = 5, grid = g1,
    seed = 21, fit_seed = 9, nthread = 1, keep_oof = TRUE, verbose = FALSE)))
  expect_identical(dim(tune$oof$pred), c(nrow(d$smp), 1L))
  ## keep_oof changes nothing else in the result
  t0 <- suppressWarnings(suppressMessages(povmap::xgb_tune(
    y ~ x1 + x2, smp_data = d$smp, smp_weights = "wt", domains = "dom", transformation = "arcsin", folds = 5, grid = g1,
    seed = 21, fit_seed = 9, nthread = 1, verbose = FALSE)))
  t1 <- tune; t1$oof <- NULL
  expect_identical(t1, t0)
  reused <- pool_fit(d, area_pool = "out_of_fold", pool_tune = tune)
  fresh <- pool_fit(d, area_pool = "out_of_fold", pool_folds = tune$fold_assignment, pool_fit_seed = 9)
  expect_true(reused$pools$area_from_tune); expect_false(fresh$pools$area_from_tune)
  expect_equal(reused$pools$configs[[1]]$resid_domains, fresh$pools$configs[[1]]$resid_domains, tolerance = 1e-10)
  expect_equal(reused$pools$configs[[1]]$wt_domains, fresh$pools$configs[[1]]$wt_domains, tolerance = 1e-12)
  ## and the sub-area pool still takes its own fits, alongside
  both <- pool_fit(d, area_pool = "out_of_fold", pool_tune = tune, sub_pool = "out_of_fold")
  expect_equal(both$pools$configs[[1]]$resid_domains, reused$pools$configs[[1]]$resid_domains, tolerance = 1e-12)
  expect_false(isTRUE(all.equal(both$pools$configs[[1]]$resid_sub_domains, reused$pools$configs[[1]]$resid_sub_domains)))
  ## mismatches stop with the reason
  expect_error(pool_fit(d, area_pool = "out_of_fold", pool_tune = tune, hp = list(nrounds = 16, max_depth = 2, eta = 0.3, min_child_weight = 1,
               lambda = 1, colsample_bytree = 0.6, colsample_bylevel = 1, colsample_bynode = 1, subsample = 0.6)), "hyperparameters")
  expect_error(pool_fit(d, area_pool = "out_of_fold", pool_tune = tune, transformation = "no"), "transformation")
  d2 <- d; d2$smp$y[3] <- min(d2$smp$y[3] + 0.01, 1)
  expect_error(pool_fit(d2, area_pool = "out_of_fold", pool_tune = tune), "outcome")
  expect_error(pool_fit(d, area_pool = "out_of_fold", pool_tune = tune, variance_y = "vy"), "variance_y")
  expect_error(pool_fit(d, pool_tune = tune), "area_pool")
  expect_error(pool_fit(d, area_pool = "out_of_fold", pool_tune = t0), "keep_oof")
  d3 <- d; d3$pop$x1 <- d3$pop$x1 + 0.1
  expect_error(pool_fit(d3, area_pool = "out_of_fold", pool_tune = tune), "differ")
})

test_that("DRC Stage 1b scratch hook reproduced by xgb(area_pool/sub_pool = out of fold) on access_water, fold 1, B = 100", {
  ## inst/scripts/oof_pool_drc_reproduction.R re-runs the held-out territories of CV fold 1 (seed 123) with the package
  ## options and compares with the saved scratch-hook job; this checks its small output file when present.
  f <- "H:/david/SAE/DRC/data/candidate_estimates/drc_retune_20261010/oof_pool_repro_access_water_f01.csv"
  skip_on_cran()
  skip_if(!file.exists(f), "reproduction output not available")
  r <- utils::read.csv(f)
  expect_equal(r$n_territories[1], 21)
  expect_lt(max(r$max_abs_diff), 1e-8)                    # Mean, Mean_bench, Lower, Upper, Lower_bench, Upper_bench
  expect_lt(max(r$max_abs_diff_vs_insample_point, na.rm = TRUE), 1e-10)   # point estimates equal the in-sample job's
})
