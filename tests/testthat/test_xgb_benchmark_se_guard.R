# Regression test for the guard on benchmark_target_se in xgb().
#
# When a benchmark group's direct estimate is on the boundary of its range (a proportion of exactly 0 or 1),
# its design-based standard error is exactly 0. With perturb_benchmark = TRUE the perturbation of the benchmark
# target is rnorm(B, 0, sd = SE), so a zero SE adds no uncertainty for that group and the intervals of its
# domains collapse towards the target. This was silent. xgb() now stops, telling the caller to supply a
# boundary-corrected SE (Agresti-Coull for a proportion at 0 or 1); the package never substitutes one, and the
# stop can be overridden deliberately with allow_zero_benchmark_se = TRUE (zero only: negative and
# non-finite values always stop).
#
# Tested here:
#   1. a zero SE stops, naming the group and Agresti-Coull, before any bootstrap work;
#   2. negative, NA, NaN and infinite SEs stop, also with the override;
#   3. perturb_benchmark = FALSE does not look at the SEs, so a zero target SE is fine;
#   4. allow_zero_benchmark_se = TRUE proceeds with a warning, and a zero-SE group's replicates are unperturbed
#      (benchmarked variance equals the unperturbed variance) while the other groups' are not;
#   5. the internal Horvitz-Thompson fallback (no benchmark_target_se) is guarded the same way: a group whose
#      sampled units all share one outcome has an internal SE of exactly 0;
#   6. valid SEs pass silently.

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

fit_guard <- function(d, se, perturb = TRUE, ...) {
  args <- list(fixed = y ~ x1 + x2, smp_data = d$smp, smp_weights = "wt", pop_data = d$pop, pop_weights = "npop",
               domains = "dom", sub_domains = "sub", transformation = "arcsin", bootstrap = TRUE, B = 30, L = 5, cpus = 1,
               seed = 3, nrounds = 15, eta = 0.3, max_depth = 2, colsample_bylevel = 1, colsample_bynode = 1, subsample = 1,
               benchmark = "Mean", benchmark_level = "reg", benchmark_type = "logit_raking",
               perturb_benchmark = perturb, benchmark_target_se = se)
  do.call(povmap::xgb, c(args, list(...)))
}
quiet <- function(expr) suppressMessages(suppressWarnings(expr))

test_that("a zero benchmark_target_se stops with advice, naming the group", {
  d <- make_guard_data()
  expect_error(quiet(fit_guard(d, c(R1 = 0.02, R2 = 0, R3 = 0.02))), "zero for 1 benchmark group")
  expect_error(quiet(fit_guard(d, c(R1 = 0.02, R2 = 0, R3 = 0.02))), "R2")
  expect_error(quiet(fit_guard(d, c(R1 = 0.02, R2 = 0, R3 = 0.02))), "Agresti-Coull")
  expect_error(quiet(fit_guard(d, c(R1 = 0.02, R2 = 0, R3 = 0.02))), "allow_zero_benchmark_se")
  expect_error(quiet(fit_guard(d, c(R1 = 0, R2 = 0, R3 = 0.02))), "zero for 2 benchmark group")
})

test_that("negative and non-finite standard errors always stop", {
  d <- make_guard_data()
  for (bad in list(c(R1 = 0.02, R2 = -0.01, R3 = 0.02), c(R1 = 0.02, R2 = NA, R3 = 0.02),
                   c(R1 = 0.02, R2 = NaN, R3 = 0.02), c(R1 = 0.02, R2 = Inf, R3 = 0.02))) {
    expect_error(quiet(fit_guard(d, bad)), "finite and non-negative")
    expect_error(quiet(fit_guard(d, bad, allow_zero_benchmark_se = TRUE)), "finite and non-negative")   # the override covers zero only
  }
})

test_that("nothing is checked when the perturbation is off", {
  d <- make_guard_data()
  fit <- quiet(fit_guard(d, c(R1 = 0.02, R2 = 0, R3 = 0.02), perturb = FALSE))
  expect_s3_class(fit, "xgb")
  expect_true(all(is.finite(fit$ind$Mean_bench)))
  fit2 <- quiet(fit_guard(d, c(R1 = 0.02, R2 = -1, R3 = NA), perturb = FALSE))           # not used, so not validated
  expect_s3_class(fit2, "xgb")
})

test_that("allow_zero_benchmark_se = TRUE proceeds deliberately, with a warning, and the zero group is unperturbed", {
  d <- make_guard_data()
  se <- c(R1 = 0.05, R2 = 0, R3 = 0.05)
  expect_warning(fit <- suppressMessages(fit_guard(d, se, allow_zero_benchmark_se = TRUE)), "R2")
  doms_r2 <- unique(d$pop$dom[d$pop$reg == "R2"])
  v <- fit$var$Var_bench[match(doms_r2, fit$var$Domain)]
  vu <- fit$Var_bench_unperturbed$Var_bench_unperturbed[match(doms_r2, fit$Var_bench_unperturbed$Domain)]
  expect_equal(v, vu, tolerance = 1e-12)                                                   # nothing was added for R2
  doms_r1 <- unique(d$pop$dom[d$pop$reg == "R1"])
  expect_true(all(fit$var$Var_bench[match(doms_r1, fit$var$Domain)] >
                  fit$Var_bench_unperturbed$Var_bench_unperturbed[match(doms_r1, fit$Var_bench_unperturbed$Domain)]))
})

test_that("the internal Horvitz-Thompson fallback is guarded the same way", {
  d <- make_guard_data()
  fit <- quiet(fit_guard(d, NULL))                                                         # ordinary data: SE > 0, passes
  expect_s3_class(fit, "xgb")
  ## every sampled unit of R2 has the same outcome: the internal HT SE of R2 is exactly 0
  d2 <- d; d2$smp$y[d2$smp$reg == "R2"] <- 0.4
  expect_error(quiet(fit_guard(d2, NULL)), "internal Horvitz-Thompson standard error")
  expect_error(quiet(fit_guard(d2, NULL)), "zero for 1 benchmark group")
  expect_error(quiet(fit_guard(d2, NULL)), "R2")
  expect_error(quiet(fit_guard(d2, NULL)), "Agresti-Coull")
  ## deliberate override: proceeds with a warning
  expect_warning(f2 <- suppressMessages(fit_guard(d2, NULL, allow_zero_benchmark_se = TRUE)), "internal Horvitz-Thompson")
  expect_s3_class(f2, "xgb")
  ## no perturbation, no check
  expect_s3_class(quiet(fit_guard(d2, NULL, perturb = FALSE)), "xgb")
  ## supplying positive SEs for every group avoids the internal estimate altogether
  expect_s3_class(quiet(fit_guard(d2, c(R1 = 0.02, R2 = 0.02, R3 = 0.02))), "xgb")
})

test_that("valid standard errors pass silently", {
  d <- make_guard_data()
  expect_no_warning(suppressMessages(fit_guard(d, c(R1 = 0.02, R2 = 0.01, R3 = 0.03))))
  ## a very small but positive SE is valid
  expect_s3_class(quiet(fit_guard(d, c(R1 = 1e-8, R2 = 0.01, R3 = 0.03))), "xgb")
})
