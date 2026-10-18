# What megb ships to parallel bootstrap workers: the benchmark framework carries only the columns benchmark_fn reads,
# and none of it changes a number.
test_that(".megb_keep_cols keeps the named columns in order and leaves other objects alone", {
  df <- data.frame(a = 1:3, b = 4:6, c = 7:9)
  expect_identical(povmap:::.megb_keep_cols(df, c("c", "a", "zz")), df[, c("a", "c")])
  expect_identical(povmap:::.megb_keep_cols(df, character(0)), df[, character(0), drop = FALSE])
  expect_identical(povmap:::.megb_keep_cols(1:3, "a"), 1:3)
  if (requireNamespace("data.table", quietly = TRUE)) {
    dt <- data.table::as.data.table(df)
    expect_identical(povmap:::.megb_keep_cols(dt, "a"), dt)
  }
})

test_that(".megb_lean_framework drops the census columns benchmark_fn never reads and nothing else", {
  pop <- data.frame(dom = rep(letters[1:4], each = 5), reg = rep(c("A", "B"), each = 10), w = 1:20 + 0.5,
                    x1 = rnorm(20), x2 = rnorm(20), x3 = rnorm(20))
  smp <- data.frame(dom = rep(letters[1:4], each = 2), reg = rep(c("A", "B"), each = 4), y = runif(8), sw = 1, x1 = rnorm(8), x2 = rnorm(8))
  fwk <- list(Y_smp = smp$y, smp_data = smp, pop_data = pop, domains = "dom", smp_domains = "dom", outcome = "y",
              smp_weights = "sw", pop_weights = "w", benchmark_weights = "sw", aggregate_to = NULL, n_pop = c(5L, 5L, 5L, 5L))
  lean <- povmap:::.megb_lean_framework(fwk, "reg")
  expect_identical(names(lean), names(fwk))
  expect_identical(names(lean$pop_data), c("dom", "reg", "w"))
  expect_identical(names(lean$smp_data), c("dom", "reg", "y", "sw"))
  for (nm in setdiff(names(fwk), c("pop_data", "smp_data"))) expect_identical(lean[[nm]], fwk[[nm]])
  expect_identical(lean$pop_data$w, fwk$pop_data$w)
  # the benchmark itself is unchanged on the lean framework
  est <- c(0.2, 0.3, 0.4, 0.5); bm <- c(A = 0.35, B = 0.45)
  for (bt in c("ratio", "logit_raking"))
    expect_identical(
      benchmark_xgb_level(list(ind = data.frame(Domain = letters[1:4], Mean = est)), lean, y ~ x1, bm, bt, "reg"),
      benchmark_xgb_level(list(ind = data.frame(Domain = letters[1:4], Mean = est)), fwk, y ~ x1, bm, bt, "reg"))
})

test_that("the lean framework is much smaller to serialise when the census has many unused columns", {
  pop <- data.frame(dom = rep(1:50, each = 40), reg = rep(1:5, each = 400), w = runif(2000))
  for (j in 1:60) pop[[paste0("x", j)]] <- rnorm(2000)
  fwk <- list(smp_data = data.frame(dom = 1:50, reg = rep(1:5, each = 10), y = runif(50), sw = 1), pop_data = pop,
              domains = "dom", smp_domains = "dom", outcome = "y", smp_weights = "sw", pop_weights = "w", benchmark_weights = "sw")
  lean <- povmap:::.megb_lean_framework(fwk, "reg")
  expect_lt(length(serialize(lean, NULL)), length(serialize(fwk, NULL)) / 10)
})

test_that("megb with cpus = 2 (full refit, benchmark and perturbation) is reproducible and shares the point estimates and perturbation draws with cpus = 1", {
  # The replicate streams differ between cpus = 1 (the session RNG) and cpus > 1 (furrr L'Ecuyer streams), as in v0.6.2, so the
  # bootstrap intervals are compared across two cpus = 2 runs, not with cpus = 1.
  skip_on_cran()
  skip_if_not_installed("furrr")
  cfg <- megb_bench_cfgs$log_ratio
  d <- megb_bench_data()
  tgt <- tapply(d$smp$y, d$smp$reg, mean); tgt <- stats::setNames(as.numeric(tgt) * c(1.05, 0.95, 1.02), names(tgt))
  run <- function(cpus) megb_bench_run(cfg, bootstrap_refit = "full", B = 4, em_maxit = 3L, cpus = cpus,
                                       perturb_benchmark = TRUE, benchmark_target_se = 0.05 * tgt)
  a <- run(1L); b <- run(2L); b2 <- run(2L)
  expect_identical(b$ind, b2$ind)
  expect_identical(b$CI, b2$CI)
  expect_identical(a$ind$Mean, b$ind$Mean)
  expect_identical(a$ind$Mean_bench, b$ind$Mean_bench)
  expect_identical(a$benchmark_perturbation, b$benchmark_perturbation)
})
