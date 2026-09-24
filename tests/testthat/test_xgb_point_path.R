# Regression tests for the point-estimate path of xgb().
#
# point_estim_xgb() forms the point estimate as the fitted prediction plus a
# sub-area and an area residual draw, averaged over L draws. The sub-area pool
# is centred by construction. The area pool was never centred: its expectation
# under its draw probability is the overall weighted mean residual, which shifted
# every domain by the same amount, and the shift did not fall with L. In the
# Colombia migration SAE study it cost the xgb arm 0.67 percent in RMSE on its
# own, with a further 0.80 percent of Monte Carlo noise at L = 100. The switch
# that centres the pools, center_residuals, reached only the bootstrap branch.
#
# Fixed behaviour, tested here:
#   1. both pools the point path draws from are mean zero under their own draw
#      probabilities;
#   2. with transformation = "no" the point estimate IS the deterministic
#      population-weighted aggregate of the fitted predictions, with no draw.

make_data <- function(seed = 7) {
  set.seed(seed)
  nd <- 30; ns <- 12
  pop <- expand.grid(sa = seq_len(ns), dom = seq_len(nd))
  pop$sub <- sprintf("%02d_%02d", pop$dom, pop$sa)
  pop$x1 <- rnorm(nrow(pop)); pop$x2 <- runif(nrow(pop))
  u <- rnorm(nd, 0, 0.4)[pop$dom]
  pop$npop <- round(runif(nrow(pop), 50, 500))
  pop$p <- plogis(-2 + 0.8 * pop$x1 + 0.6 * pop$x2 + u)
  pop$dom <- sprintf("%02d", pop$dom)
  ## sample: 18 domains, 6 sub-areas each, unequal weights
  sd_ <- sort(sample(unique(pop$dom), 18))
  smp <- do.call(rbind, lapply(sd_, function(d) {
    s <- pop[pop$dom == d, ][sample(ns, 6), ]
    s$y <- rbinom(nrow(s), 25, s$p) / 25
    s$wt <- exp(rnorm(nrow(s), 0, 0.7)) * 100
    s }))
  list(pop = pop[, c("sub", "dom", "x1", "x2", "npop")],
       smp = smp[, c("sub", "dom", "x1", "x2", "y", "wt")])
}

capture_point_path <- function(expr) {
  env <- new.env()
  suppressMessages(trace("point_estim_xgb", where = asNamespace("povmap"), print = FALSE,
                         exit = bquote(assign("out", returnValue(), envir = .(env)))))
  on.exit(suppressMessages(untrace("point_estim_xgb", where = asNamespace("povmap"))))
  res <- force(expr)
  list(res = res, pe = env$out)
}

fit_xgb <- function(d, transformation, L = 20, seed = 1) {
  povmap::xgb(fixed = y ~ x1 + x2, smp_data = d$smp, smp_weights = "wt",
              pop_data = d$pop, pop_weights = "npop", domains = "dom",
              sub_domains = "sub", transformation = transformation,
              bootstrap = FALSE, L = L, seed = seed, nrounds = 15, eta = 0.3,
              max_depth = 2, colsample_bylevel = 1, colsample_bynode = 1)
}

test_that("both residual pools the point path draws from are mean zero", {
  skip_if_not_installed("povmap")
  skip_if_not_installed("xgboost")
  d <- make_data()
  for (tr in c("arcsin", "no")) {
    cp <- capture_point_path(fit_xgb(d, tr))
    pe <- cp$pe
    # sub-area pool, drawn with prob = the raw survey weight
    expect_equal(stats::weighted.mean(pe$resid_sub_domains, pe$wt_sub_domains), 0,
                 tolerance = 1e-12, info = tr)
    # area pool used by the point path, drawn with prob = the summed domain weight
    expect_equal(stats::weighted.mean(pe$resid_domains_point, pe$wt_domains), 0,
                 tolerance = 1e-12, info = tr)
    # and the uncentred pool is NOT mean zero on this fit, so the test bites
    expect_gt(abs(stats::weighted.mean(pe$resid_domains, pe$wt_domains)), 1e-4)
  }
})

test_that("the point-path area draw uses the centred pool", {
  skip_if_not_installed("povmap")
  b <- paste(deparse(povmap:::point_estim_xgb), collapse = "\n")
  expect_true(grepl("resid_domains_point[sample(", b, fixed = TRUE))
  expect_false(grepl("area_draw = resid_domains[sample(", b, fixed = TRUE))
})

test_that("untransformed point estimate equals the deterministic aggregate", {
  skip_if_not_installed("povmap")
  skip_if_not_installed("xgboost")
  d <- make_data()
  r <- fit_xgb(d, "no", L = 20)
  hat <- stats::predict(r$model, as.matrix(d$pop[, c("x1", "x2")]))
  det <- tapply(hat * d$pop$npop, d$pop$dom, sum) / tapply(d$pop$npop, d$pop$dom, sum)
  got <- setNames(r$ind$Mean, r$ind$Domain)
  expect_equal(unname(got[names(det)]), unname(as.numeric(det)), tolerance = 1e-12)
  expect_equal(r$ind$Mean_agg, r$ind$Mean)
  # no draw: the estimate cannot depend on L
  r2 <- fit_xgb(d, "no", L = 3)
  expect_identical(r2$ind$Mean, r$ind$Mean)
})
