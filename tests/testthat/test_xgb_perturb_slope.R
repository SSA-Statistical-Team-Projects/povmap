# Regression test for the delta-method slope of the benchmark-target perturbation (found by the DRC audit).
#
# xgb() turns the rate-scale SE of a benchmark group's target into an increment on the transformed scale with
# dz/dp evaluated at the group's mean benchmarked rate. For arcsin, dz/dp = 1 / (2 sqrt(p (1 - p))) is unbounded at 0 and 1.
# A group benchmarked to a target of exactly 0 had its slope taken at about 1e-9, i.e. about 100 (a finite difference over 1e-4).
# With a zero SE the product was 0, which hid it; with a positive boundary-corrected SE (Agresti-Coull) the perturbation became
# a fraction of a radian to a radian (100 x 0.005 = 0.5): half the
# replicates ended at 0, the rest spread over (0, 1), and the upper bounds of the group's domains reached 0.1 to 0.6.
# The slope is now taken at no less than the target SE from either bound (the rate is clamped into [se, 1 - se], arcsin only).
#
# Tested here:
#   1. .perturb_slope(): identical to the old computation wherever the clamp does not bind (mean rate at least one SE from
#      the boundary), clamped where it does; se = 0 / NA leave it alone; log and sqrt are untouched;
#   2. through xgb(), with keep_replicates: the realised slope (perturbed minus unperturbed benchmarked replicate on the
#      arcsin scale, over the noise draw) equals the unclamped slope for groups away from the boundary, and the clamped
#      one for a group pinned at zero with a positive SE;
#   3. that pinned group's domains get sensible widths (upper bounds of the order of the target SE, not 0.1 to 0.6).

old_slope <- function(pbar, tf, hi) {                     # the computation before the clamp
  lo <- 1e-9; pbar <- min(max(pbar, lo), hi); eps <- 1e-4
  p_hi <- min(pbar + eps, hi); p_lo <- max(pbar - eps, lo)
  (tf(p_hi)$y - tf(p_lo)$y) / (p_hi - p_lo)
}
asin_tf <- function(y) list(y = asin(sqrt(y)))

test_that(".perturb_slope is unchanged where the clamp does not bind", {
  se <- 0.01
  for (p in c(0.011, 0.02, 0.1, 0.5, 0.9, 0.98, 0.989))
    expect_identical(povmap:::.perturb_slope(p, se, "arcsin", asin_tf), old_slope(p, asin_tf, 1 - 1e-9), info = p)
  ## exactly one SE from the boundary: still the old value
  expect_identical(povmap:::.perturb_slope(se, se, "arcsin", asin_tf), old_slope(se, asin_tf, 1 - 1e-9))
  expect_identical(povmap:::.perturb_slope(1 - se, se, "arcsin", asin_tf), old_slope(1 - se, asin_tf, 1 - 1e-9))
})

test_that(".perturb_slope clamps the rate to [se, 1 - se] for arcsin", {
  se <- 0.005
  s_at_se <- 1 / (2 * sqrt(se * (1 - se)))
  for (p in c(0, 1e-9, 0.001, 0.004999)) {
    expect_equal(povmap:::.perturb_slope(p, se, "arcsin", asin_tf), s_at_se, tolerance = 1e-3, info = p)
    expect_gte(old_slope(p, asin_tf, 1 - 1e-9), 0.999 * povmap:::.perturb_slope(p, se, "arcsin", asin_tf))   # the old slope was at least as large
  }
  expect_gt(old_slope(0, asin_tf, 1 - 1e-9), 90)                    # the old slope at a target of 0 (a finite difference over 1e-4: about 100)
  expect_lt(povmap:::.perturb_slope(0, se, "arcsin", asin_tf), 10)  # now about 7
  for (p in c(1, 1 - 1e-9, 0.9995)) expect_equal(povmap:::.perturb_slope(p, se, "arcsin", asin_tf), s_at_se, tolerance = 1e-3, info = p)
  ## a huge SE: the slope is taken at 0.5
  expect_equal(povmap:::.perturb_slope(0, 0.7, "arcsin", asin_tf), 1, tolerance = 1e-6)
  ## no SE to go on: unchanged
  expect_identical(povmap:::.perturb_slope(0, 0, "arcsin", asin_tf), old_slope(0, asin_tf, 1 - 1e-9))
  expect_identical(povmap:::.perturb_slope(0, NA_real_, "arcsin", asin_tf), old_slope(0, asin_tf, 1 - 1e-9))
})

test_that(".perturb_slope leaves the other transformations alone", {
  lg <- function(y) list(y = log(y)); sq <- function(y) list(y = sqrt(y))
  for (p in c(0.001, 0.05, 0.2, 3)) for (se in c(0, 0.01, 0.3)) {
    expect_identical(povmap:::.perturb_slope(p, se, "log", lg), old_slope(p, lg, Inf), info = paste(p, se))
    expect_identical(povmap:::.perturb_slope(p, se, "sqrt", sq), old_slope(p, sq, Inf), info = paste(p, se))
  }
})

make_slope_data <- function(seed = 51) {
  set.seed(seed)
  nd <- 24; ns <- 8
  pop <- expand.grid(sa = seq_len(ns), dom = seq_len(nd))
  pop$sub <- sprintf("%02d_%02d", pop$dom, pop$sa)
  pop$reg <- sprintf("R%d", (pop$dom - 1) %/% 8 + 1)
  pop$x1 <- rnorm(nrow(pop)); pop$x2 <- runif(nrow(pop))
  u <- rnorm(nd, 0, 0.4)[pop$dom]
  pop$npop <- round(runif(nrow(pop), 50, 500))
  pop$p <- plogis(ifelse(pop$reg == "R3", -2.2, 0) + 0.6 * pop$x1 + 0.4 * pop$x2 + u)
  pop$dom <- sprintf("%02d", pop$dom)
  pop$r2 <- as.numeric(pop$reg == "R2"); pop$r3 <- as.numeric(pop$reg == "R3")      # region dummies, as the state dummies in production
  sd_ <- unlist(lapply(split(unique(pop$dom), (seq_len(nd) - 1) %/% 8), function(d) sort(sample(d, 5))))
  smp <- do.call(rbind, lapply(sd_, function(d) pop[pop$dom == d, ][sample(ns, 4), ]))
  smp$y <- rbinom(nrow(smp), 25, smp$p) / 25
  smp$y[smp$reg == "R2"] <- 0                                       # R2: every sampled outcome 0, so the target is exactly 0
  smp$wt <- exp(rnorm(nrow(smp), 0, 0.4)) * 100
  list(pop = pop[, c("sub", "dom", "reg", "x1", "x2", "r2", "r3", "npop")], smp = smp[, c("sub", "dom", "reg", "x1", "x2", "r2", "r3", "y", "wt")])
}
fit_slope <- function(d, se = c(R1 = 0.02, R2 = 0.005, R3 = 0.01)) {
  suppressWarnings(suppressMessages(povmap::xgb(
    fixed = y ~ x1 + x2 + r2 + r3, smp_data = d$smp, smp_weights = "wt", pop_data = d$pop, pop_weights = "npop", domains = "dom", sub_domains = "sub",
    transformation = "arcsin", bootstrap = TRUE, B = 150, L = 5, cpus = 1, seed = 3, nrounds = 20, eta = 0.3, max_depth = 2,
    colsample_bylevel = 1, colsample_bynode = 1, subsample = 1, benchmark = "Mean", benchmark_level = "reg", benchmark_type = "logit_raking",
    perturb_benchmark = TRUE, benchmark_target_se = se, keep_replicates = TRUE)))
}

test_that("through xgb(): realised slopes, and sensible widths for a group pinned at zero", {
  d <- make_slope_data()
  se <- c(R1 = 0.02, R2 = 0.005, R3 = 0.01)
  fit <- fit_slope(d, se)
  r <- fit$replicates
  to <- function(y) asin(sqrt(pmin(pmax(y, 1e-9), 1 - 1e-9)))
  realised <- function(g) {                                           # perturbed minus unperturbed, on the arcsin scale, per unit noise
    doms <- unique(d$pop$dom[d$pop$reg == g]); j <- match(doms, r$domains)
    zu <- to(r$benchmarked_unperturbed[, j, drop = FALSE]); zp <- to(r$benchmarked[, j, drop = FALSE])
    noise <- r$benchmark_perturbation[, g]
    ok <- abs(noise) > 1e-4 & zp > 1e-6 & zp < pi / 2 - 1e-6 & zu > 1e-6 & zu < pi / 2 - 1e-6   # not clipped by the back-transform
    list(slope = median(((zp - zu) / noise)[ok & matrix(TRUE, nrow(zp), ncol(zp))]),
         pbar = mean(r$benchmarked_unperturbed[, j]), n = sum(ok))
  }
  for (g in c("R1", "R3")) {                                          # away from the boundary: the clamp does not bind
    x <- realised(g)
    expect_gt(x$pbar, 3 * se[[g]]); expect_lt(x$pbar, 1 - 3 * se[[g]])
    expect_gt(x$n, 100)
    expect_equal(x$slope, old_slope(x$pbar, asin_tf, 1 - 1e-9), tolerance = 1e-4, info = g)
  }
  x2 <- realised("R2")                                                # pinned at zero: slope at the SE, not at 1e-9
  expect_lt(x2$pbar, 0.01)
  expect_gt(x2$n, 20)
  expect_equal(x2$slope, 1 / (2 * sqrt(se[["R2"]] * (1 - se[["R2"]]))), tolerance = 0.02)
  expect_lt(x2$slope, 10)                                             # the old slope here was about 16,000
  ## widths: the target SE is 0.005, so upper bounds of the order of 0.005 to 0.02, not 0.1 to 0.6
  doms2 <- unique(d$pop$dom[d$pop$reg == "R2"]); k <- match(doms2, fit$CI$Domain)
  expect_lt(max(fit$CI$Upper_bench[k]), 0.05)
  expect_lt(max(fit$replicates$benchmarked[, match(doms2, r$domains)]), 0.1)      # no replicate is thrown across (0, 1)
  expect_true(all(fit$CI$Lower_bench[k] >= 0 & fit$CI$Upper_bench[k] > fit$CI$Lower_bench[k]))
})
