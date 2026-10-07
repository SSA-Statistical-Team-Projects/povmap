# External benchmark targets by benchmark_level group.
#
# A survey's published province direct estimates are computed on the full sample,
# in each EA's design province. The model sample can differ (EAs without a
# location, EAs whose grid cell falls in another province), so benchmarking to the
# sample mean of the modelled EAs ("Mean") does not reproduce the published
# figures. benchmark_xgb_level() benchmarks to a numeric vector named by the
# benchmark_level groups, but xgb_check1() accepted a numeric benchmark only as a
# national Mean/Head_Count of length 1 or 2, so that path was unreachable.
#
# Tested here: with benchmark = a named vector of group targets, the
# population-weighted aggregate of the benchmarked domain estimates equals each
# group's target; malformed target vectors stop with a clear message.

make_level_data <- function(seed = 3) {
  set.seed(seed)
  nd <- 30; ns <- 12
  pop <- expand.grid(sa = seq_len(ns), dom = seq_len(nd))
  pop$sub <- sprintf("%02d_%02d", pop$dom, pop$sa)
  pop$reg <- sprintf("R%d", (pop$dom - 1) %/% 10 + 1)
  pop$x1 <- rnorm(nrow(pop)); pop$x2 <- runif(nrow(pop))
  u <- rnorm(nd, 0, 0.4)[pop$dom]
  pop$npop <- round(runif(nrow(pop), 50, 500))
  pop$p <- plogis(-1 + 0.8 * pop$x1 + 0.6 * pop$x2 + u)
  pop$dom <- sprintf("%02d", pop$dom)
  sd_ <- unlist(lapply(split(unique(pop$dom), (seq_len(nd) - 1) %/% 10),
                       function(d) sort(sample(d, 6))))
  smp <- do.call(rbind, lapply(sd_, function(d) {
    s <- pop[pop$dom == d, ][sample(ns, 6), ]
    s$y <- rbinom(nrow(s), 25, s$p) / 25
    s$wt <- exp(rnorm(nrow(s), 0, 0.7)) * 100
    s }))
  list(pop = pop[, c("sub", "dom", "reg", "x1", "x2", "npop")],
       smp = smp[, c("sub", "dom", "reg", "x1", "x2", "y", "wt")])
}

fit_level <- function(d, benchmark, type = "ratio") {
  suppressWarnings(povmap::xgb(
    fixed = y ~ x1 + x2, smp_data = d$smp, smp_weights = "wt",
    pop_data = d$pop, pop_weights = "npop", domains = "dom", sub_domains = "sub",
    transformation = "no", bootstrap = FALSE, L = 5, seed = 1, nrounds = 15,
    eta = 0.3, max_depth = 2, colsample_bylevel = 1, colsample_bynode = 1,
    subsample = 1, benchmark = benchmark, benchmark_level = "reg",
    benchmark_type = type))
}

group_aggregate <- function(d, est) {
  dp <- tapply(d$pop$npop, d$pop$dom, sum)
  rg <- tapply(d$pop$reg, d$pop$dom, `[`, 1)
  e  <- setNames(est$ind$Mean_bench, as.character(est$ind$Domain))[names(dp)]
  sapply(split(names(dp), rg[names(dp)]), function(ds) sum(e[ds] * dp[ds]) / sum(dp[ds]))
}

test_that("benchmarked domain estimates aggregate to each group's external target", {
  skip_if_not_installed("xgboost")
  d <- make_level_data()
  targets <- c(R1 = 0.31, R2 = 0.47, R3 = 0.58)
  for (type in c("ratio", "logit_raking")) {
    est <- fit_level(d, targets, type)
    agg <- group_aggregate(d, est)
    expect_equal(unname(agg[names(targets)]), unname(targets), tolerance = 1e-8, info = type)
  }
  # and these are not the sample means, so the test bites
  smp_mean <- sapply(split(d$smp, d$smp$reg), function(s) stats::weighted.mean(s$y, s$wt))
  expect_gt(max(abs(smp_mean[names(targets)] - targets)), 0.02)
})

test_that("malformed group targets stop with a clear message", {
  skip_if_not_installed("xgboost")
  d <- make_level_data()
  expect_error(fit_level(d, c(0.3, 0.4, 0.5)), "named by")
  expect_error(fit_level(d, c(R1 = 0.3, R2 = 0.4)), "no target for the reg group\\(s\\) R3")
  expect_error(fit_level(d, c(R1 = 0.3, R2 = NA, R3 = 0.5)), "R2 are NA")
  expect_error(fit_level(d, c(R1 = 0.3, R1 = 0.4, R3 = 0.5)), "unique name")
})

test_that("a national Mean benchmark is still accepted without benchmark_level", {
  skip_if_not_installed("xgboost")
  d <- make_level_data()
  est <- suppressWarnings(povmap::xgb(
    fixed = y ~ x1 + x2, smp_data = d$smp, smp_weights = "wt",
    pop_data = d$pop, pop_weights = "npop", domains = "dom", sub_domains = "sub",
    transformation = "no", bootstrap = FALSE, L = 5, seed = 1, nrounds = 15,
    eta = 0.3, max_depth = 2, colsample_bylevel = 1, colsample_bynode = 1,
    subsample = 1, benchmark = c(Mean = 0.45), benchmark_type = "ratio"))
  expect_true("Mean_bench" %in% names(est$ind))
})
