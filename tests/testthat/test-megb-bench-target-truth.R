# megb(bench_target = "truth"): each bootstrap replicate is benchmarked to T* = the simulated truth of its benchmark
# group (aggregated with the weights the benchmark aggregates the estimate with) plus a draw from N(0, SE_g^2), so the
# provincial error of the benchmarked estimate is that draw alone. The point estimate is unchanged.

grp_dom <- c(d1 = "A", d2 = "A", d3 = "B", d4 = "B", d5 = "B")
w_dom   <- c(d1 = 10, d2 = 30, d3 = 5, d4 = 5, d5 = 10)

test_that(".megb_target_star aggregates the domain truths with the domain totals and adds the draw (G x B)", {
  set.seed(3)
  tau <- matrix(runif(5 * 4, 0.2, 0.6), 5, 4)
  draws <- matrix(c(0.01, -0.02, 0.03, 0, 0.1, 0.2, -0.1, 0.05), 4, 2, dimnames = list(NULL, c("A", "B")))
  ts <- povmap:::.megb_target_star(tau, names(grp_dom), grp_dom, w_dom, draws)
  expect_equal(dim(ts), c(2L, 4L)); expect_identical(rownames(ts), c("A", "B"))
  tA <- colSums(tau[1:2, ] * w_dom[1:2]) / sum(w_dom[1:2]); tB <- colSums(tau[3:5, ] * w_dom[3:5]) / sum(w_dom[3:5])
  expect_equal(unname(ts["A", ]), tA + draws[, "A"]); expect_equal(unname(ts["B", ]), tB + draws[, "B"])
  # domains given in another order than the group / weight vectors
  ts2 <- povmap:::.megb_target_star(tau[5:1, ], rev(names(grp_dom)), grp_dom, w_dom, draws)
  expect_equal(ts2, ts)
})

test_that(".megb_target_star clamps the target (arcsin: [1e-9, 1 - 1e-9]; no: unclamped); SE 0 returns the truth", {
  tau <- matrix(c(0.97, 0.98, 0.02, 0.01, 0.03), 5, 3)
  draws <- matrix(c(0.2, 0.2, 0.2, -0.2, -0.2, -0.2), 3, 2, dimnames = list(NULL, c("A", "B")))
  b <- povmap:::.megb_bench_bounds("arcsin")
  ts <- povmap:::.megb_target_star(tau, names(grp_dom), grp_dom, w_dom, draws, b[["lo"]], b[["hi"]])
  expect_true(all(ts >= 1e-9 & ts <= 1 - 1e-9)); expect_equal(unname(ts["A", ]), rep(1 - 1e-9, 3)); expect_equal(unname(ts["B", ]), rep(1e-9, 3))
  expect_equal(unname(povmap:::.megb_bench_bounds("log")), c(1e-9, Inf)); expect_equal(unname(povmap:::.megb_bench_bounds("no")), c(-Inf, Inf))
  zero <- matrix(0, 3, 2, dimnames = list(NULL, c("A", "B")))
  tz <- povmap:::.megb_target_star(tau, names(grp_dom), grp_dom, w_dom, zero)
  expect_equal(unname(tz["A", ]), unname(colSums(tau[1:2, ] * w_dom[1:2]) / sum(w_dom[1:2])))
})

test_that("the target-error draws are the perturbation's draws (same stream), and the caller's RNG state is left alone", {
  se <- c(A = 0.03, B = 0.05)
  set.seed(42); s0 <- .Random.seed
  a <- povmap:::.megb_bench_draws(se, c("A", "B"), 50, 9, FALSE, "bench_target = \"truth\"")$draws
  expect_identical(.Random.seed, s0)
  e <- list(tau_b_bench = matrix(0.5, 6, 50), tau_star_orig_bench = matrix(0.5, 6, 50), bench_keep = NULL, MSE_bench_estimates = NULL)
  p <- povmap:::.megb_perturb_bench(e, se, FALSE, B = 50, seed = 9, dom_group = c("A", "A", "A", "B", "B", "B"), transformation = "no",
                                   transform_outcome = no_transform, back_transform_outcome = no_transform_back)
  expect_identical(a, p$benchmark_perturbation)
})

test_that("bench_target = \"truth\" needs what it needs and excludes perturb_benchmark", {
  cfg <- megb_bench_cfgs$no_ratio; se <- c(A = 0.1, B = 0.1, C = 0.1)
  expect_error(megb_bench_run(cfg, bench_target = "truth", benchmark_target_se = se, perturb_benchmark = TRUE), "replaces perturb_benchmark")
  expect_error(megb_bench_run(cfg, bench_target = "truth"), "needs benchmark_target_se")
  expect_error(megb_bench_run(cfg, bench_target = "truth", benchmark_target_se = se, benchmark_level = NULL, benchmark = NULL), "needs a benchmark, a benchmark_level")
  expect_error(megb_bench_run(cfg, bench_target = "truth", benchmark_target_se = se, pop_weights = NULL), "needs pop_weights")
  expect_error(megb_bench_run(cfg, bench_target = "truth", benchmark_target_se = c(A = 0.1, B = 0.1)), "missing 1 of 3")
  expect_error(megb_bench_run(cfg, bench_target = "truth", benchmark_target_se = c(A = 0.1, B = 0, C = 0.1)), "zero")
  expect_error(megb_bench_run(cfg, bench_target = "other"), "should be one of")
})

# Run megb with the bootstrap collation captured, to read the benchmarked replicates and the truth they are compared with.
run_captured <- function(cfg, ...) {
  env <- new.env(); orig <- povmap:::.mse_megb_collate
  out <- testthat::with_mocked_bindings(megb_bench_run(cfg, ...),
                                        .mse_megb_collate = function(...) { r <- orig(...); assign("res", r, envir = env); r }, .package = "povmap")
  list(out = out, mse = get("res", envir = env))
}

# Group error of the benchmarked replicates: sum_d w_d (est_d - truth_d) / sum_d w_d, a G x B' matrix.
group_error <- function(cap, d = megb_bench_data()) {
  m <- cap$mse; doms <- as.character(m$MSE_estimates[[1]])
  w <- as.numeric(tapply(d$pop$w, as.character(d$pop$dom), sum)[doms])
  g <- d$pop$reg[match(doms, as.character(d$pop$dom))]
  err <- m$tau_b_bench - m$tau_star_orig_bench
  t(vapply(sort(unique(g)), function(k) colSums(err[g == k, , drop = FALSE] * w[g == k]) / sum(w[g == k]), numeric(ncol(err))))
}

for (nm in names(megb_bench_cfgs)) {
  cfg <- megb_bench_cfgs[[nm]]
  test_that(paste0("truth target (", nm, "): the group error of the benchmarked replicates is the drawn target error"), {
    se <- c(A = 0.02, B = 0.03, C = 0.04)
    cap <- run_captured(cfg, bench_target = "truth", benchmark_target_se = se)
    ge <- group_error(cap); keep <- cap$mse$bench_keep; if (is.null(keep)) keep <- seq_len(ncol(ge))
    draws <- cap$out$benchmark_target_draws
    expect_equal(dim(draws), c(6L, 3L)); expect_identical(colnames(draws), c("A", "B", "C")); expect_identical(cap$out$bench_target, "truth")
    for (g in c("A", "B", "C")) expect_equal(unname(ge[g, ]), unname(draws[keep, g]), tolerance = 1e-6)
  })

  test_that(paste0("truth target (", nm, "): SE = 0 gives zero group error; the point estimate is unchanged"), {
    zero <- c(A = 0, B = 0, C = 0)
    cap <- run_captured(cfg, bench_target = "truth", benchmark_target_se = zero, allow_zero_benchmark_se = TRUE)
    expect_lt(max(abs(group_error(cap))), 1e-6)
    base <- megb_bench_run(cfg)
    expect_identical(cap$out$ind$Mean, base$ind$Mean); expect_identical(cap$out$ind$Mean_bench, base$ind$Mean_bench)
  })
}

test_that("truth target: under 'fixed' the group error varies with the simulated truth, under 'truth' its spread is the target SE only", {
  cfg <- megb_bench_cfgs$no_ratio; se <- c(A = 0.02, B = 0.02, C = 0.02)
  fx <- group_error(run_captured(cfg, bench_target = "fixed")); tr <- group_error(run_captured(cfg, bench_target = "truth", benchmark_target_se = se))
  expect_true(all(apply(fx, 1, sd) > 0))
  # the fixed-target group error is T - truth_g: its spread is the simulated truth's; with a small SE the truth-target error is far smaller
  expect_lt(mean(apply(tr, 1, sd)), mean(apply(fx, 1, sd)))
})

test_that("truth target: the arcsin target stays in [0, 1] and the replicates are raked inside the range, large SE included", {
  cfg <- megb_bench_cfgs$arcsin_logit
  cap <- run_captured(cfg, bench_target = "truth", benchmark_target_se = c(A = 0.6, B = 0.6, C = 0.6))
  expect_true(all(cap$mse$tau_b_bench >= 0 & cap$mse$tau_b_bench <= 1))
  expect_true(all(is.finite(cap$out$CI_bench$Lower)) || all(is.finite(cap$out$CI$Lower_bench)))
})

test_that("truth target works with the leaves_only and full bootstrap refits and gives the same draws and points", {
  cfg <- megb_bench_cfgs$no_ratio; se <- c(A = 0.05, B = 0.05, C = 0.05)
  for (rf in c("leaves_only", "full")) {
    r <- tryCatch(megb_bench_run(cfg, bootstrap_refit = rf, bench_target = "truth", benchmark_target_se = se, B = 4,
                                 gradient_params = list(eta = 0.3, max_depth = 2, subsample = 1, nrounds = 20), em_maxit = 2),
                  error = function(e) e)
    expect_false(inherits(r, "error"), info = if (inherits(r, "error")) paste(rf, conditionMessage(r)))
    if (!inherits(r, "error")) {
      expect_equal(dim(r$benchmark_target_draws), c(4L, 3L)); expect_true(all(is.finite(r$CI$Lower_bench)))
      expect_identical(r$ind$Mean_bench, megb_bench_run(cfg, bootstrap_refit = rf, B = 4, em_maxit = 2,
                                                         gradient_params = list(eta = 0.3, max_depth = 2, subsample = 1, nrounds = 20))$ind$Mean_bench)
    }
  }
})

test_that("truth target: full refit with cpus = 2 has the same draws and points as cpus = 1", {
  skip_if(Sys.getenv("POVMAP_TEST_CPUS") != "1", "set POVMAP_TEST_CPUS=1 to run (starts 2 worker processes)")
  skip_on_cran()
  cfg <- megb_bench_cfgs$no_ratio; se <- c(A = 0.05, B = 0.05, C = 0.05)
  gp <- list(eta = 0.3, max_depth = 2, subsample = 1, nrounds = 20)
  a <- megb_bench_run(cfg, bootstrap_refit = "full", bench_target = "truth", benchmark_target_se = se, B = 4, em_maxit = 2, gradient_params = gp, cpus = 1)
  b <- megb_bench_run(cfg, bootstrap_refit = "full", bench_target = "truth", benchmark_target_se = se, B = 4, em_maxit = 2, gradient_params = gp, cpus = 2)
  expect_identical(a$benchmark_target_draws, b$benchmark_target_draws); expect_identical(a$ind$Mean_bench, b$ind$Mean_bench)
})

test_that("the default (no bench_target) is bench_target = \"fixed\", and the result then has no truth-target elements", {
  cfg <- megb_bench_cfgs$log_ratio
  a <- megb_bench_run(cfg); b <- megb_bench_run(cfg, bench_target = "fixed")
  expect_identical(a$ind, b$ind); expect_identical(a$CI, b$CI); expect_identical(a$var, b$var)
  expect_null(a$benchmark_target_draws); expect_null(a$bench_target)
})

# ---- arcsin_backtransform ----------------------------------------------------------------------------------------
test_that(".megb_make_corrected_bt: 'plain' is sin^2 whatever the variance; 'exact' and other transformations are unchanged", {
  th <- c(-0.2, 0.1, 0.6, 1.2, 1.7)
  pl <- povmap:::.megb_make_corrected_bt(0.09, "arcsin", arcsin_plain = TRUE)
  ex <- povmap:::.megb_make_corrected_bt(0.09, "arcsin")
  expect_identical(pl(th, 0.3), arcsin_transform_back(th)); expect_identical(pl(th), arcsin_transform_back(th))
  expect_false(attr(pl, "uses_re_var")); expect_true(attr(ex, "uses_re_var"))
  expect_identical(attr(ex, "uses_re_var"), TRUE)
  expect_equal(ex(th, 0.2), povmap:::arcsin_expected_back(th, 0.09 + 0.2))
  expect_false(isTRUE(all.equal(ex(th, 0.2), pl(th, 0.2))))
  expect_identical(povmap:::.megb_re_var(pl, NULL, c("a", "b")), c(0, 0))
  for (tr in c("log", "poisson", "no")) {   # the flag is for arcsin only
    a <- povmap:::.megb_make_corrected_bt(0.09, tr, arcsin_plain = TRUE); b <- povmap:::.megb_make_corrected_bt(0.09, tr)
    expect_identical(a(th, 0.1), b(th, 0.1)); expect_identical(attr(a, "uses_re_var"), attr(b, "uses_re_var"))
  }
})

test_that("arcsin_backtransform = 'exact' is the default; 'plain' equals REStat's in-process override of the factory", {
  cfg <- megb_bench_cfgs$arcsin_logit
  a <- megb_bench_run(cfg); b <- megb_bench_run(cfg, arcsin_backtransform = "exact")
  expect_identical(a$ind, b$ind); expect_identical(a$CI, b$CI); expect_identical(a$var, b$var); expect_null(a$arcsin_backtransform)
  p <- megb_bench_run(cfg, arcsin_backtransform = "plain")
  expect_identical(p$arcsin_backtransform, "plain"); expect_false(isTRUE(all.equal(a$ind$Mean, p$ind$Mean)))
  # REStat's scripts/_ct4_megb_task_bt.R, arm P: the factory is replaced by one returning arcsin_transform_back(hat_t)
  orig <- povmap:::.megb_make_corrected_bt
  restat <- function(sigma2_e, transformation, ls_lambda = NULL, arcsin_plain = FALSE) {
    f0 <- orig(sigma2_e, transformation, ls_lambda); if (!identical(transformation, "arcsin")) return(f0)
    g <- function(hat_t, re_var = 0) arcsin_transform_back(hat_t); attr(g, "uses_re_var") <- FALSE
    environment(g) <- list2env(list(), parent = asNamespace("povmap")); g }
  r <- testthat::with_mocked_bindings(megb_bench_run(cfg), .megb_make_corrected_bt = restat, .package = "povmap")
  expect_identical(p$ind, r$ind); expect_identical(p$CI, r$CI); expect_identical(p$var, r$var)
  expect_identical(p$ind_bench, r$ind_bench); expect_identical(p$CI_bench, r$CI_bench)
})

test_that("arcsin_backtransform = 'plain' with benchmarking and the truth target, and it is inert for other transformations", {
  cfg <- megb_bench_cfgs$arcsin_logit; se <- c(A = 0.02, B = 0.02, C = 0.02)
  r <- megb_bench_run(cfg, arcsin_backtransform = "plain", bench_target = "truth", benchmark_target_se = se)
  expect_true(all(is.finite(r$CI$Lower_bench)) && all(r$ind$Mean_bench >= 0 & r$ind$Mean_bench <= 1))
  l <- megb_bench_cfgs$log_ratio
  expect_identical(megb_bench_run(l)$ind, megb_bench_run(l, arcsin_backtransform = "plain")$ind)
  expect_null(megb_bench_run(l, arcsin_backtransform = "plain")$arcsin_backtransform)
  expect_error(megb_bench_run(cfg, arcsin_backtransform = "series"), "should be one of")
})
