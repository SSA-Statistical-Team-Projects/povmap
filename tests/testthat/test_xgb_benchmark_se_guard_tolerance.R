# Follow-up to test_xgb_benchmark_se_guard.R: what that file leaves open.
#
#   1. Alignment. The SEs are matched to the benchmark groups by NAME. A reordered vector gives the same fit as an
#      ordered one, the guard names the right group, and an unnamed vector is refused. (Without name matching, a
#      reordered SE vector would pass the check and perturb the wrong province.)
#   2. Tolerance. The threshold at which an SE counts as zero scales with the outcome: a floating-point near-zero is
#      about 1e-16 times the outcome, so an absolute 1e-10 missed it for large outcomes (an internal Horvitz-Thompson
#      SE of identical units was 5e-10 for an outcome of 1e7), and a supplied 1e-17 passed. It must still never catch
#      a legitimate SE.
#   3. No bootstrap, no perturbation: the SEs are not used, so they are not validated.

make_guard_data <- function(seed = 11) {
  set.seed(seed)
  nd <- 24; ns <- 8
  pop <- expand.grid(sa = seq_len(ns), dom = seq_len(nd))
  pop$sub <- sprintf("%02d_%02d", pop$dom, pop$sa)
  pop$reg <- sprintf("R%d", (pop$dom - 1) %/% 8 + 1)
  pop$x1 <- rnorm(nrow(pop)); pop$x2 <- runif(nrow(pop))
  u <- rnorm(nd, 0, 0.4)[pop$dom]
  pop$npop <- round(runif(nrow(pop), 50, 500))
  pop$p <- plogis(-0.3 + 0.8 * pop$x1 + 0.6 * pop$x2 + u)
  pop$dom <- sprintf("%02d", pop$dom)
  sd_ <- unlist(lapply(split(unique(pop$dom), (seq_len(nd) - 1) %/% 8), function(d) sort(sample(d, 5))))
  smp <- do.call(rbind, lapply(sd_, function(d) pop[pop$dom == d, ][sample(ns, 4), ]))
  smp$y <- rbinom(nrow(smp), 25, smp$p) / 25
  smp$wt <- exp(rnorm(nrow(smp), 0, 0.5)) * 100
  list(pop = pop[, c("sub", "dom", "reg", "x1", "x2", "npop")], smp = smp[, c("sub", "dom", "reg", "x1", "x2", "y", "wt")])
}
fit_guard <- function(d, se, perturb = TRUE, bootstrap = TRUE, transformation = "arcsin", benchmark_type = "logit_raking", ...) {
  args <- list(fixed = y ~ x1 + x2, smp_data = d$smp, smp_weights = "wt", pop_data = d$pop, pop_weights = "npop",
               domains = "dom", sub_domains = "sub", transformation = transformation, bootstrap = bootstrap, B = 30, L = 5, cpus = 1,
               seed = 3, nrounds = 15, eta = 0.3, max_depth = 2, colsample_bylevel = 1, colsample_bynode = 1, subsample = 1,
               benchmark = "Mean", benchmark_level = "reg", benchmark_type = benchmark_type,
               perturb_benchmark = perturb, benchmark_target_se = se)
  do.call(povmap::xgb, c(args, list(...)))
}
quiet <- function(expr) suppressMessages(suppressWarnings(expr))

test_that("SEs are matched to benchmark groups by name: a reordered vector cannot pass the check and perturb the wrong group", {
  d <- make_guard_data()
  ## the zero is still found, and attributed to R2, whatever the order
  expect_error(quiet(fit_guard(d, c(R3 = 0.03, R1 = 0.02, R2 = 0))), "zero for 1 benchmark group.*: R2")
  expect_error(quiet(fit_guard(d, c(R2 = 0, R3 = 0.03, R1 = 0.02))), "zero for 1 benchmark group.*: R2")
  ## an ordered and a reordered vector give the same fit ...
  f1 <- quiet(fit_guard(d, c(R1 = 0.05, R2 = 0.001, R3 = 0.04)))
  f2 <- quiet(fit_guard(d, c(R3 = 0.04, R2 = 0.001, R1 = 0.05)))
  expect_equal(f1$CI, f2$CI)
  expect_equal(f1$var, f2$var)
  ## ... and the comparison has power: attaching the same SEs to different groups changes it
  f3 <- quiet(fit_guard(d, c(R1 = 0.001, R2 = 0.05, R3 = 0.04)))
  expect_false(isTRUE(all.equal(f1$var, f3$var)))
  ## a vector with no names cannot be matched, so it is refused
  expect_error(quiet(fit_guard(d, c(0.02, 0.02, 0.02))), "named vector")
})

test_that("the zero tolerance scales with the outcome and never catches a legitimate SE", {
  chk <- povmap:::.check_benchmark_se
  tolf <- povmap:::.benchmark_se_tol
  ## a floating-point near-zero is caught at any scale of the outcome ...
  se <- c(A = 4e-10, B = 0.5); yb <- c(A = 1e7, B = 0.5)         # 4e-10 is above an absolute 1e-10
  expect_error(chk(se, "x", FALSE, tol = tolf(se, yb)), "zero for 1 benchmark group.*: A")
  se <- c(A = 5e-8, B = 2e9); yb <- c(A = 1e9, B = 3e9)
  expect_error(chk(se, "x", FALSE, tol = tolf(se, yb)), "zero for 1 benchmark group.*: A")
  ## ... including a supplied 1e-17 beside ordinary SEs, for a group whose own mean is 0
  se <- c(A = 1e-17, B = 0.02, C = 0.03); yb <- c(A = 0, B = 0.2, C = 0.3)
  expect_error(chk(se, "x", FALSE, tol = tolf(se, yb)), "zero for 1 benchmark group.*: A")
  ## exact zeros and all-zero outcomes (mean 0) are still caught
  expect_error(chk(c(A = 0, B = 0.02), "x", FALSE, tol = tolf(c(A = 0, B = 0.02), c(A = 0, B = 0.2))), "zero for 1 benchmark group.*: A")
  ## a legitimate SE is never caught: tiny outcome scales, and SEs several orders of magnitude below the others
  se <- c(A = 2e-11, B = 4e-11); yb <- c(A = 5e-9, B = 6e-9)     # an outcome measured in tiny units
  expect_silent(chk(se, "x", FALSE, tol = tolf(se, yb)))
  se <- c(A = 1e-5, B = 0.05); yb <- c(A = 0.3, B = 0.3)
  expect_silent(chk(se, "x", FALSE, tol = tolf(se, yb)))
  se <- c(A = 8.9e4, B = 3e5); yb <- c(A = 1.8e6, B = 1.9e6)     # DRC pcexp scale
  expect_silent(chk(se, "x", FALSE, tol = tolf(se, yb)))
  ## a tolerance vector is matched to the SEs by name, not by position
  expect_error(chk(c(A = 1e-3, B = 0.5), "x", FALSE, tol = c(B = 1e-12, A = 1e-2)), "zero for 1 benchmark group.*: A")
})

test_that("end to end: a supplied near-zero SE and a large-scale identical-outcome group are stopped", {
  d <- make_guard_data()
  expect_error(quiet(fit_guard(d, c(R1 = 0.02, R2 = 1e-17, R3 = 0.03))), "zero for 1 benchmark group.*: R2")
  expect_s3_class(quiet(fit_guard(d, c(R1 = 0.02, R2 = 1e-5, R3 = 0.03))), "xgb")       # small but legitimate: passes
  ## an outcome of the order of 1e9 (a monetary amount in small units): the group of identical outcomes has a
  ## floating-point internal SE that an absolute 1e-10 would not catch
  d2 <- d
  d2$smp$y <- d2$smp$y * 1e9
  d2$smp$y[d2$smp$reg == "R2"] <- 4e8
  r2 <- d2$smp$reg == "R2"
  se_r2 <- NA
  for (s in 1:300) {
    set.seed(s); d2$smp$wt[r2] <- exp(rnorm(sum(r2), 0, 0.5)) * 100
    se_r2 <- sqrt(povmap:::ht_var_weighted_mean(d2$smp$y[r2], d2$smp$wt[r2], rep("R2", sum(r2))))
    if (se_r2 > 1e-10) break
  }
  expect_gt(se_r2, 1e-10)                                                               # beyond the old absolute tolerance
  expect_error(quiet(fit_guard(d2, NULL, transformation = "no", benchmark_type = "ratio")),
               "internal Horvitz-Thompson standard error.*zero for 1 benchmark group.*: R2")
})

test_that("with bootstrap = FALSE the SEs are not used, so they are not checked", {
  d <- make_guard_data()
  fit <- quiet(fit_guard(d, c(R1 = 0.02, R2 = 0, R3 = 0.02), bootstrap = FALSE))
  expect_s3_class(fit, "xgb")
  d2 <- d; d2$smp$y[d2$smp$reg == "R2"] <- 0.4
  expect_s3_class(quiet(fit_guard(d2, NULL, bootstrap = FALSE)), "xgb")
  ## the same data with the bootstrap on still stops
  expect_error(quiet(fit_guard(d2, NULL)), "zero for 1 benchmark group.*: R2")
})
