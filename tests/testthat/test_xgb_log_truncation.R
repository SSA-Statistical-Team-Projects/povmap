# A log-transformed outcome (an inequality index such as theil) is positive, but on the natural interval scale the
# bootstrap quantiles can put a lower bound below zero. xgb() truncated proportions to [0, 1] and poisson / gamma
# outcomes at 0, but not log-transformed ones. The lower bounds (Lower and Lower_bench) are now truncated at 0 for
# transformation = "log" as well.

make_log_data <- function(seed = 5) {
  set.seed(seed)
  nd <- 24; ns <- 10
  pop <- expand.grid(sa = seq_len(ns), dom = seq_len(nd))
  pop$sub <- sprintf("%02d_%02d", pop$dom, pop$sa)
  pop$reg <- sprintf("R%d", (pop$dom - 1) %/% 8 + 1)
  pop$x1 <- rnorm(nrow(pop)); pop$x2 <- runif(nrow(pop))
  u <- rnorm(nd, 0, 0.5)[pop$dom]
  pop$npop <- round(runif(nrow(pop), 50, 500))
  pop$p <- plogis(-2.5 + 0.9 * pop$x1 + 0.5 * pop$x2 + u)        # a low index: many wards near 0
  pop$dom <- sprintf("%02d", pop$dom)
  sd_ <- unlist(lapply(split(unique(pop$dom), (seq_len(nd) - 1) %/% 8), function(d) sort(sample(d, 5))))
  smp <- do.call(rbind, lapply(sd_, function(d) pop[pop$dom == d, ][sample(ns, 4), ]))
  set.seed(9); smp$y <- pmax(0.01, 0.08 * rlnorm(nrow(smp), 0, 1.2))   # a skewed, index-like outcome with a long right tail
  smp$wt <- exp(rnorm(nrow(smp), 0, 0.5)) * 100
  list(pop = pop[, c("sub", "dom", "reg", "x1", "x2", "npop")], smp = smp[, c("sub", "dom", "reg", "x1", "x2", "y", "wt")])
}
fit_log <- function(d, ...) suppressWarnings(suppressMessages(povmap::xgb(
  fixed = y ~ x1 + x2, smp_data = d$smp, smp_weights = "wt", pop_data = d$pop, pop_weights = "npop", domains = "dom", sub_domains = "sub",
  transformation = "log", bootstrap = TRUE, B = 40, L = 5, cpus = 1, seed = 3, nrounds = 20, eta = 0.3, max_depth = 2,
  colsample_bylevel = 1, colsample_bynode = 1, subsample = 1, benchmark = "Mean", benchmark_level = "reg", benchmark_type = "ratio", ...)))

test_that("log: natural-scale lower bounds are truncated at 0", {
  d <- make_log_data()
  for (m in c("basic", "plus")) {
    f <- fit_log(d, interval_method = m, interval_scale = "natural", keep_replicates = FALSE)
    expect_true(all(f$CI$Lower >= 0), info = m)
    expect_true(all(f$CI$Lower_bench >= 0), info = m)
    ## the truncation is active on this data (before the change, with "basic", 7 wards had a negative Lower_bench and 1 a negative Lower)
    if (m == "basic") expect_true(any(f$CI$Lower_bench == 0))
    expect_true(all(f$CI$Lower_bench[f$CI$Lower_bench > 0] < f$CI$Upper_bench[f$CI$Lower_bench > 0]), info = m)
    expect_true(all(f$CI$Lower_bench <= f$CI$Upper_bench), info = m)
  }
})

test_that("log: the transformed interval scale is unchanged (positive by construction)", {
  d <- make_log_data()
  f <- fit_log(d, interval_scale = "transformed")
  expect_true(all(f$CI$Lower > 0) && all(f$CI$Lower_bench > 0))
})

test_that("log: xgb() stops on an outcome that is not strictly positive, so truncating at 0 is always valid", {
  d <- make_log_data()
  for (v in c(0, -0.3)) {
    d2 <- d; d2$smp$y[1] <- v
    expect_error(fit_log(d2), "strictly greater than 0")
  }
})
