# megb(perturb_benchmark = TRUE): the benchmark-target perturbation of xgb, applied to the benchmarked
# bootstrap replicates after the bootstrap, on the transformed scale, truth side unperturbed.

fake_mse_est <- function(D = 6, Bk = 40, keep = NULL, seed = 1) {
  set.seed(seed)
  M <- matrix(runif(D * Bk, 0.2, 0.7), D, Bk)
  Tr <- matrix(runif(D * Bk, 0.2, 0.7), D, Bk)
  list(tau_b_bench = M, tau_star_orig_bench = Tr, bench_keep = keep,
       MSE_bench_estimates = data.frame(Domain = sprintf("d%d", seq_len(D)), Mean = rowMeans((Tr - M)^2)))
}
grp6 <- c("A", "A", "A", "B", "B", "B")
se2 <- c(A = 0.03, B = 0.05)

test_that(".megb_perturb_bench adds the pre-drawn noise group by group (transformation 'no': exactly)", {
  e <- fake_mse_est()
  p <- povmap:::.megb_perturb_bench(e, se2, FALSE, B = 40, seed = 9, dom_group = grp6, transformation = "no",
                                   transform_outcome = no_transform, back_transform_outcome = no_transform_back)
  expect_identical(p$tau_b_bench_unperturbed, e$tau_b_bench)
  expect_identical(p$tau_star_orig_bench, e$tau_star_orig_bench)               # truth side untouched
  expect_equal(dim(p$benchmark_perturbation), c(40L, 2L)); expect_identical(colnames(p$benchmark_perturbation), c("A", "B"))
  for (g in c("A", "B")) {
    r <- which(grp6 == g)
    expect_equal(p$tau_b_bench[r, ] - e$tau_b_bench[r, ], matrix(p$benchmark_perturbation[, g], length(r), 40, byrow = TRUE))
  }
  expect_equal(p$MSE_bench_estimates$Mean, rowMeans((e$tau_star_orig_bench - p$tau_b_bench)^2))
})

test_that(".megb_perturb_bench on the transformed scale: slope times draw, back-transformed, stays in range (arcsin)", {
  e <- fake_mse_est(); e$tau_b_bench[1:3, ] <- 0.999; e$tau_b_bench[4:6, ] <- 0.001
  p <- povmap:::.megb_perturb_bench(e, c(A = 0.08, B = 0.08), FALSE, B = 40, seed = 3, dom_group = grp6, transformation = "arcsin",
                                   transform_outcome = arcsin_transform, back_transform_outcome = arcsin_transform_back)
  expect_true(all(p$tau_b_bench >= 0 & p$tau_b_bench <= 1))                      # no clamp needed
  r <- which(grp6 == "A")
  dz <- povmap:::.perturb_slope(mean(e$tau_b_bench[r, ]), 0.08, "arcsin", arcsin_transform)
  expected <- arcsin_transform_back(arcsin_transform(pmin(pmax(as.vector(e$tau_b_bench[r, ]), 1e-9), 1 - 1e-9))$y +
                                    as.vector(matrix(p$benchmark_perturbation[, "A"], length(r), 40, byrow = TRUE)) * dz)
  expect_equal(as.vector(p$tau_b_bench[r, ]), as.vector(expected))
})

test_that(".megb_perturb_bench uses the draws of the kept replicates only", {
  keep <- c(2L, 5L, 7L, 8L, 30L)
  e <- fake_mse_est(Bk = length(keep), keep = keep)
  p <- povmap:::.megb_perturb_bench(e, se2, FALSE, B = 40, seed = 9, dom_group = grp6, transformation = "no",
                                   transform_outcome = no_transform, back_transform_outcome = no_transform_back)
  expect_equal(p$tau_b_bench[1, ] - e$tau_b_bench[1, ], p$benchmark_perturbation[keep, "A"])
})

test_that("the draws depend only on the seed, have the target SD, and leave the caller's RNG state alone", {
  e <- fake_mse_est(Bk = 4)
  args <- list(mse_est = e, benchmark_target_se = se2, allow_zero = FALSE, dom_group = grp6, transformation = "no",
               transform_outcome = no_transform, back_transform_outcome = no_transform_back)
  set.seed(42); s0 <- .Random.seed
  a <- do.call(povmap:::.megb_perturb_bench, c(args, list(B = 20000, seed = 11)))
  expect_identical(.Random.seed, s0)
  b <- do.call(povmap:::.megb_perturb_bench, c(args, list(B = 20000, seed = 11)))
  d <- do.call(povmap:::.megb_perturb_bench, c(args, list(B = 20000, seed = 12)))
  expect_identical(a$benchmark_perturbation, b$benchmark_perturbation)
  expect_false(identical(a$benchmark_perturbation, d$benchmark_perturbation))
  expect_equal(apply(a$benchmark_perturbation, 2, sd), se2, tolerance = 0.03)
  expect_equal(unname(colMeans(a$benchmark_perturbation)), c(0, 0), tolerance = 0.01, scale = 1)
})

test_that("benchmark_target_se is validated like xgb: names, missing groups, bad and zero values", {
  e <- fake_mse_est(Bk = 4)
  f <- function(se, allow = FALSE) povmap:::.megb_perturb_bench(e, se, allow, B = 4, seed = 1, dom_group = grp6, transformation = "no",
                                                                 transform_outcome = no_transform, back_transform_outcome = no_transform_back)
  expect_error(f(c(0.03, 0.05)), "named vector")
  expect_error(f(c(A = 0.03)), "missing 1 of 2")
  expect_error(f(c(A = 0.03, B = -0.1)), "finite and non-negative")
  expect_error(f(c(A = 0.03, B = NA_real_)), "finite and non-negative")
  expect_error(f(c(A = 0.03, B = 0)), "is zero for 1 benchmark group")
  expect_warning(z <- f(c(A = 0.03, B = 0), allow = TRUE), "allow_zero_benchmark_se = TRUE")
  expect_true(all(z$benchmark_perturbation[, "B"] == 0))
  expect_equal(z$tau_b_bench[4:6, ], e$tau_b_bench[4:6, ])                      # a zero SE adds nothing there
})

test_that("perturb_benchmark = FALSE reproduces the v0.6.1 megb output exactly (benchmarked, all transformations tested)", {
  fx <- testthat::test_path("fixtures", "megb_bench_identity_v061.rds")
  skip_if_not(file.exists(fx), "fixture missing")
  ref <- readRDS(fx)
  for (nm in names(megb_bench_cfgs)) for (mt in c("var", "mse")) {
    m <- megb_bench_run(megb_bench_cfgs[[nm]], mse_type = mt)
    r <- ref[[paste(nm, mt, sep = "|")]]
    expect_identical(m$ind, r$ind); expect_identical(m$var, r$var); expect_identical(m$CI, r$CI)
    expect_null(m$Var_bench_unperturbed); expect_null(m$benchmark_perturbation)
    m2 <- megb_bench_run(megb_bench_cfgs[[nm]], mse_type = mt, perturb_benchmark = FALSE, benchmark_target_se = c(A = 0.1, B = 0.1, C = 0.1))
    expect_identical(m2$CI, r$CI)
  }
})

test_that("megb(perturb_benchmark = TRUE) leaves the point estimates and the bootstrap unchanged and widens the benchmarked intervals", {
  se3 <- c(A = 0.06, B = 0.06, C = 0.06)
  for (nm in names(megb_bench_cfgs)) {
    cfg <- megb_bench_cfgs[[nm]]
    se_use <- if (cfg$tr == "no") se3 * 5 else se3
    m0 <- megb_bench_run(cfg, B = 40)
    m1 <- megb_bench_run(cfg, B = 40, perturb_benchmark = TRUE, benchmark_target_se = se_use)
    expect_identical(m1$ind, m0$ind)                                              # point estimates
    expect_identical(m1$CI[, c("Domain", "Lower", "Upper")], m0$CI[, c("Domain", "Lower", "Upper")])   # unbenchmarked interval
    expect_identical(m1$var$Mean, m0$var$Mean)
    u <- m1$Var_bench_unperturbed
    expect_identical(sort(u$Domain), sort(as.character(m0$var$Domain)))
    expect_equal(u$Var_bench_unperturbed[match(m0$var$Domain, u$Domain)], m0$var$Var_bench)   # same bootstrap
    expect_gt(mean(m1$var$Var_bench), mean(m0$var$Var_bench))
    w1 <- m1$CI$Upper_bench - m1$CI$Lower_bench; w0 <- m0$CI$Upper_bench - m0$CI$Lower_bench
    expect_gt(mean(w1), mean(w0))
    expect_equal(dim(m1$benchmark_perturbation), c(40L, 3L))
  }
})

test_that("megb(perturb_benchmark = TRUE) also runs with mse_type = 'mse'", {
  se3 <- c(A = 0.06, B = 0.06, C = 0.06)
  m0 <- megb_bench_run(megb_bench_cfgs$arcsin_logit, mse_type = "mse", B = 20)
  m1 <- megb_bench_run(megb_bench_cfgs$arcsin_logit, mse_type = "mse", B = 20, perturb_benchmark = TRUE, benchmark_target_se = se3)
  expect_gt(mean(m1$var$Var_bench), mean(m0$var$Var_bench))
  expect_equal(m1$Var_bench_unperturbed$Var_bench_unperturbed[match(m0$var$Domain, m1$Var_bench_unperturbed$Domain)], m0$var$Var_bench)
})

test_that("a zero perturbation (allow_zero_benchmark_se) reproduces the unperturbed result", {
  z <- c(A = 0, B = 0, C = 0)
  m0 <- megb_bench_run(megb_bench_cfgs$log_ratio, B = 20)
  m1 <- suppressWarnings(megb_bench_run(megb_bench_cfgs$log_ratio, B = 20, perturb_benchmark = TRUE, benchmark_target_se = z,
                                        allow_zero_benchmark_se = TRUE))
  expect_equal(m1$CI, m0$CI)
  expect_equal(m1$var, m0$var)
})

test_that("perturb_benchmark = TRUE needs its inputs", {
  cfg <- megb_bench_cfgs$arcsin_logit
  se3 <- c(A = 0.06, B = 0.06, C = 0.06)
  expect_error(megb_bench_run(cfg, perturb_benchmark = TRUE), "benchmark_target_se")
  expect_error(megb_bench_run(cfg, perturb_benchmark = TRUE, benchmark_target_se = se3, mse = FALSE), "mse = TRUE")
  expect_error(megb_bench_run(cfg, perturb_benchmark = TRUE, benchmark_target_se = se3, benchmark_level = NULL), "benchmark_level")
  expect_error(megb_bench_run(cfg, perturb_benchmark = NA), "TRUE or FALSE")
  expect_error(megb_bench_run(cfg, perturb_benchmark = TRUE, benchmark_target_se = c(A = 0.06, B = 0.06)), "missing 1 of 3")
})

test_that("the lower clamp applies to arcsin, log and poisson only; log.shift is not clamped", {
  lam <- 2
  tf <- function(y) list(y = log(y + lam), shift = NULL); bf <- function(y) exp(y) - lam
  e <- fake_mse_est(); e$tau_b_bench[1:3, ] <- matrix(runif(3 * 40, -1.2, -0.3), 3, 40); e$tau_b_bench[4:6, ] <- matrix(runif(3 * 40, 0.3, 1.2), 3, 40)
  p <- povmap:::.megb_perturb_bench(e, c(A = 0.05, B = 0.05), FALSE, B = 40, seed = 3, dom_group = grp6, transformation = "log.shift",
                                   transform_outcome = tf, back_transform_outcome = bf)
  for (g in c("A", "B")) {
    r <- which(grp6 == g)
    dz <- povmap:::.perturb_slope(mean(e$tau_b_bench[r, ]), 0.05, "log.shift", tf)
    expected <- bf(tf(as.vector(e$tau_b_bench[r, ]))$y + as.vector(matrix(p$benchmark_perturbation[, g], length(r), 40, byrow = TRUE)) * dz)
    expect_equal(as.vector(p$tau_b_bench[r, ]), expected)
  }
  expect_lt(min(p$tau_b_bench[1:3, ]), 0)                                      # negative replicates are kept, not clipped to 1e-9
  e0 <- fake_mse_est(); e0$tau_b_bench[1, ] <- 0
  for (tr in c("log", "poisson")) {
    f <- if (tr == "log") list(log_transform, log_transform_back) else list(poisson_transform, poisson_transform_back)
    q <- povmap:::.megb_perturb_bench(e0, c(A = 0.05, B = 0.05), FALSE, B = 40, seed = 3, dom_group = grp6, transformation = tr,
                                     transform_outcome = f[[1]], back_transform_outcome = f[[2]])
    expect_true(all(is.finite(q$tau_b_bench[1, ])))                              # a zero is clamped to 1e-9 before the transformation
    expect_true(all(q$tau_b_bench[1, ] >= 0))
  }
})
