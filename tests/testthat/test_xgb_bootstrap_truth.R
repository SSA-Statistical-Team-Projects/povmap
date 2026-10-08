# Regression tests for the simulated truth in the benchmarked xgb bootstrap.
#
# Each replicate simulates a population: every cell gets a sub-area residual
# and its domain's area draw. The truth of a domain is the population-weighted
# mean of its simulated cells, and the interval is built from the spread of
# (benchmarked estimate - truth) across replicates.
#
# From 7f9deb0 (2026-02-11) the truth was computed as
#   collapse:::fmean(B_sub$sim_plus_area, g = B_sub$domains, w = B_sub$wts)
# but B_sub (sub_predictions) names its columns after fwk$domains and
# fwk$pop_weights, so both arguments were NULL. fmean then returned ONE
# unweighted mean of every cell, recycled across all domains. No error was
# raised. The intervals lost each domain's own area effect: for unsampled
# domains, whose refit never sees their draw, the area-level variance was
# missing altogether (DRC lfp: median width +9%, up to +54%, once corrected).
# The same NULL pattern sat in the hat_t aggregation and the case bootstrap.
#
# Tested here:
#   1. the aggregation helper refuses NULL or mis-sized grouping and weights;
#   2. in every replicate the truth differs across domains, equals each
#      domain's population-weighted mean of the simulated cells, and its
#      population-weighted mean equals that of the simulated population;
#   3. no aggregation in xgb() refers to B_sub$domains or B_sub$wts;
#   4. the case bootstrap returns one value per domain.

make_bench_data <- function(seed = 7) {
  set.seed(seed)
  nd <- 30; ns <- 12
  pop <- expand.grid(sa = seq_len(ns), dom = seq_len(nd))
  pop$sub <- sprintf("%02d_%02d", pop$dom, pop$sa)
  pop$reg <- sprintf("R%d", (pop$dom - 1) %/% 10 + 1)   # 3 benchmark regions
  pop$x1 <- rnorm(nrow(pop)); pop$x2 <- runif(nrow(pop))
  u <- rnorm(nd, 0, 0.4)[pop$dom]
  pop$npop <- round(runif(nrow(pop), 50, 500))
  pop$p <- plogis(-1 + 0.8 * pop$x1 + 0.6 * pop$x2 + u)
  pop$dom <- sprintf("%02d", pop$dom)
  ## sample: 18 of 30 domains (6 per region), 6 sub-areas each, unequal weights
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

# Record every call to the aggregation helper: its inputs and its result.
capture_boot_agg <- function(expr) {
  env <- new.env(); env$calls <- list()
  suppressMessages(trace(".boot_domain_wmean", where = asNamespace("povmap"), print = FALSE,
    exit = bquote(assign("calls", c(get("calls", envir = .(env)),
                         list(list(what = what, x = x, g = g, w = w, out = returnValue()))),
                         envir = .(env)))))
  on.exit(suppressMessages(untrace(".boot_domain_wmean", where = asNamespace("povmap"))))
  res <- force(expr)
  list(res = res, calls = env$calls)
}

fit_bench <- function(d, B = 4, bootstrap_type = "residual", benchmark = "Mean") {
  suppressWarnings(povmap::xgb(
    fixed = y ~ x1 + x2, smp_data = d$smp, smp_weights = "wt",
    pop_data = d$pop, pop_weights = "npop", domains = "dom", sub_domains = "sub",
    transformation = "no", bootstrap = TRUE, bootstrap_type = bootstrap_type,
    B = B, L = 5, cpus = 1, seed = 1, nrounds = 15, eta = 0.3, max_depth = 2,
    colsample_bylevel = 1, colsample_bynode = 1, subsample = 1,
    benchmark = benchmark, benchmark_level = if (is.null(benchmark)) NULL else "reg",
    benchmark_type = "ratio"))
}

test_that("the aggregation helper refuses NULL or mis-sized grouping and weights", {
  f <- povmap:::.boot_domain_wmean
  x <- c(1, 2, 3, 4); g <- c("a", "a", "b", "b"); w <- c(1, 3, 1, 1)
  expect_equal(unname(f(x, g, w, "t")), c(1.75, 3.5))
  expect_equal(names(f(x, g, w, "t")), c("a", "b"))
  expect_error(f(x, NULL, w, "t"), "domain grouping is NULL")
  expect_error(f(x, g, NULL, "t"), "population weight is NULL")
  expect_error(f(x, g[-1], w, "t"), "grouping has length 3, the values 4")
  expect_error(f(x, g, w[-1], "t"), "weight has length 3, the values 4")
  expect_error(f(NULL, g, w, "t"), "NULL or empty")
  # the failure the guard exists for: a $ lookup of an absent column is NULL,
  # and fmean then silently returns a single national mean
  df <- data.frame(dom = g, npop = w, v = x)
  expect_null(df$domains)
  expect_length(collapse::fmean(df$v, g = df$domains, w = df$wts), 1L)
  expect_error(f(df$v, df$domains, df$wts, "t"), "domain grouping is NULL")
})

test_that("the bootstrap truth is a per-domain, population-weighted mean of the simulated population", {
  skip_if_not_installed("xgboost")
  d <- make_bench_data()
  cp <- capture_boot_agg(fit_bench(d, B = 4))
  truth <- Filter(function(cl) cl$what == "sim_truth", cp$calls)
  expect_length(truth, 4L)                              # one per replicate
  doms <- sort(unique(d$pop$dom))
  for (b in seq_along(truth)) {
    tc <- truth[[b]]
    # grouped by the framework's domain column and weighted by its population column
    expect_equal(length(tc$g), nrow(d$pop), info = b)
    expect_setequal(unique(as.character(tc$g)), doms)
    expect_equal(sort(tc$w), sort(as.numeric(d$pop$npop)), info = b)
    # one value per domain, and they differ across domains within the replicate
    expect_equal(names(tc$out), doms, info = b)
    expect_gt(stats::sd(tc$out), 0)
    # each domain's truth is the population-weighted mean of its own simulated cells
    by_dom <- vapply(doms, function(dd) stats::weighted.mean(tc$x[tc$g == dd], tc$w[tc$g == dd]), 0)
    expect_equal(unname(tc$out), unname(by_dom), tolerance = 1e-12, info = b)
    # and the population-weighted mean of the domain truths is the simulated population mean
    pop_d <- vapply(doms, function(dd) sum(tc$w[tc$g == dd]), 0)
    expect_equal(sum(tc$out * pop_d) / sum(pop_d), stats::weighted.mean(tc$x, tc$w), tolerance = 1e-12, info = b)
  }
  # the benchmarked interval came out of a run in which all of the above held
  expect_equal(nrow(cp$res$CI), length(doms))
  expect_true(all(cp$res$CI$Upper_bench > cp$res$CI$Lower_bench))
})

test_that("no bootstrap aggregation in xgb() looks up B_sub$domains or B_sub$wts", {
  b <- paste(deparse(povmap:::xgb), collapse = "\n")
  expect_false(grepl("B_sub$domains", b, fixed = TRUE))
  expect_false(grepl("B_sub$wts", b, fixed = TRUE))
  expect_false(grepl("collapse:::fmean(B_sub", b, fixed = TRUE))
  expect_false(grepl("collapse:::fmean(x = B_sub", b, fixed = TRUE))
})

test_that("the case bootstrap aggregates each replicate to one value per domain", {
  skip_if_not_installed("xgboost")
  d <- make_bench_data()
  cp <- capture_boot_agg(fit_bench(d, B = 3, bootstrap_type = "case", benchmark = NULL))
  cs <- Filter(function(cl) cl$what == "sim (case bootstrap)", cp$calls)
  expect_length(cs, 3L)
  for (cl in cs) {
    expect_equal(names(cl$out), sort(unique(d$pop$dom)))
    expect_gt(stats::sd(cl$out), 0)
  }
})
