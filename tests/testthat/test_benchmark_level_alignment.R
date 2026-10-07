# Regression tests for domain alignment in level benchmarking.
#
# benchmark_xgb_level() bound the domain population totals to the point
# estimates by POSITION. The totals come from collapse::fsum(), whose group order
# follows collapse's `sort` option: sorted when TRUE (the package default, and
# the state of every PSOCK worker), first-appearance when FALSE (which xgb() sets
# on the master). The point estimates were likewise assumed to be in
# unique(pop_data[[domains]]) order. So whenever pop_data was not sorted by
# domain, every PARALLEL bootstrap replicate was benchmarked with domains paired
# to other domains' population weights, and its benchmarked aggregate missed its
# target. Point estimates (benchmarked on the master) were unaffected; Var_bench
# and the benchmarked intervals were not. Present in every build with
# benchmark_xgb_level(); DRC and Nigeria census tables are not sorted by domain.
#
# Tested here:
#   1. add_benchmark() hits each group's target under both sort settings, with
#      the estimates labelled in any order;
#   2. benchmark_xgb_level() matches population totals by name;
#   3. xgb() gives identical benchmarked intervals with 1 and 2 workers on
#      unsorted population data.

tiny_fwk <- function() {
  # domains deliberately NOT in sorted order in pop_data
  pop <- data.frame(dom = rep(c("D3", "D1", "D2"), each = 2), reg = "R1",
                    npop = rep(c(1000, 10, 10), each = 2))
  list(pop_data = pop, domains = "dom", pop_weights = "npop")
}

test_that("add_benchmark hits the target under both collapse sort settings", {
  old <- collapse::get_collapse("sort"); on.exit(collapse::set_collapse(sort = old))
  fwk <- tiny_fwk()
  w <- c(D3 = 2000, D1 = 20, D2 = 20)
  for (srt in c(FALSE, TRUE)) {
    collapse::set_collapse(sort = srt)
    # labelled, in pop order
    r1 <- povmap:::add_benchmark(c(0.2, 0.6, 0.6), benchmark_level = "reg", fwk = fwk,
                                 fixed = y ~ x, benchmark = c(R1 = 0.3), benchmark_type = "ratio",
                                 domains = c("D3", "D1", "D2"))
    expect_equal(sum(r1 * w[c("D3", "D1", "D2")]) / sum(w), 0.3, tolerance = 1e-12, info = srt)
    # labelled, in sorted order (as a worker's predictions come out)
    r2 <- povmap:::add_benchmark(c(0.6, 0.6, 0.2), benchmark_level = "reg", fwk = fwk,
                                 fixed = y ~ x, benchmark = c(R1 = 0.3), benchmark_type = "ratio",
                                 domains = c("D1", "D2", "D3"))
    expect_equal(sum(r2 * w[c("D1", "D2", "D3")]) / sum(w), 0.3, tolerance = 1e-12, info = srt)
    expect_equal(unname(r2), unname(r1[c(2, 3, 1)]), tolerance = 1e-12, info = srt)
    # unlabelled estimates in pop_data order still work
    r3 <- povmap:::add_benchmark(c(0.2, 0.6, 0.6), benchmark_level = "reg", fwk = fwk,
                                 fixed = y ~ x, benchmark = c(R1 = 0.3), benchmark_type = "ratio")
    expect_equal(unname(r3), unname(r1), tolerance = 1e-12, info = srt)
  }
})

test_that("labels that do not cover the population domains stop", {
  fwk <- tiny_fwk()
  expect_error(povmap:::add_benchmark(c(0.2, 0.6, 0.6), benchmark_level = "reg", fwk = fwk,
                                      fixed = y ~ x, benchmark = c(R1 = 0.3), benchmark_type = "ratio",
                                      domains = c("D3", "D1", "D9")), "cover the domains")
})

test_that("benchmark_xgb_level matches population totals to domains by name", {
  old <- collapse::get_collapse("sort"); on.exit(collapse::set_collapse(sort = old))
  collapse::set_collapse(sort = TRUE)
  b <- povmap:::benchmark_xgb_level(point_estim = list(ind = data.frame(Mean = c(0.2, 0.6, 0.6))),
                                    framework = tiny_fwk(), fixed = y ~ x, benchmark = c(R1 = 0.3),
                                    benchmark_type = "ratio", benchmark_level = "reg")
  expect_equal(sum(b$Mean_bench * c(2000, 20, 20)) / 2040, 0.3, tolerance = 1e-12)
})

test_that("benchmarked intervals do not depend on the number of workers", {
  skip_if_not_installed("xgboost")
  skip_on_cran()
  old_libs <- Sys.getenv("R_LIBS")
  Sys.setenv(R_LIBS = paste(.libPaths(), collapse = .Platform$path.sep))   # workers load this povmap
  on.exit(Sys.setenv(R_LIBS = old_libs))
  set.seed(7); nd <- 30; ns <- 12
  pop <- expand.grid(sa = seq_len(ns), dom = seq_len(nd))
  pop$sub <- sprintf("%02d_%02d", pop$dom, pop$sa)
  pop$reg <- sprintf("R%d", (pop$dom - 1) %/% 10 + 1)
  pop$x1 <- rnorm(nrow(pop)); pop$x2 <- runif(nrow(pop))
  u <- rnorm(nd, 0, 0.4)[pop$dom]
  pop$npop <- round(runif(nrow(pop), 50, 500))
  pop$p <- plogis(-1 + 0.8 * pop$x1 + 0.6 * pop$x2 + u)
  pop$dom <- sprintf("D%02d", pop$dom)
  sd_ <- unlist(lapply(split(unique(pop$dom), (seq_len(nd) - 1) %/% 10), function(d) sort(sample(d, 6))))
  smp <- do.call(rbind, lapply(sd_, function(d) {
    s <- pop[pop$dom == d, ][sample(ns, 6), ]
    s$y <- rbinom(nrow(s), 25, s$p) / 25
    s$wt <- exp(rnorm(nrow(s), 0, 0.7)) * 100
    s }))
  pop <- pop[order(-as.integer(substr(pop$dom, 2, 3)), pop$sa), ]   # domains in reverse order
  fit <- function(cpus) suppressWarnings(suppressMessages(povmap::xgb(
    fixed = y ~ x1 + x2, smp_data = smp, smp_weights = "wt",
    pop_data = pop[, c("sub", "dom", "reg", "x1", "x2", "npop")], pop_weights = "npop",
    domains = "dom", sub_domains = "sub", transformation = "no", bootstrap = TRUE,
    B = 6, L = 5, seed = 1, nrounds = 15, eta = 0.3, max_depth = 2,
    colsample_bylevel = 1, colsample_bynode = 1, subsample = 1,
    benchmark = "Mean", benchmark_level = "reg", benchmark_type = "ratio", cpus = cpus)))
  a <- fit(1); b <- fit(2)
  expect_identical(a$ind, b$ind)
  expect_equal(a$CI, b$CI, tolerance = 1e-12)
  expect_equal(a$var, b$var, tolerance = 1e-12)
})
