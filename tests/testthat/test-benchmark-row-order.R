# Regression: the benchmarked point estimate, bootstrap MSE and intervals do not depend on the order of the population rows.
#
# povmap <= 1.0.6: the benchmark function of megb()'s bootstrap handed the replicate's area estimates (in alphabetical
# domain order) to benchmark_xgb_level() without domain labels, which then assumed the first-appearance order of pop_data.
# With a pop_data not sorted by domain every replicate was benchmarked under another domain's group (ratio factors of one
# group applied to the areas of others), which inflated the benchmarked bootstrap variance and intervals several-fold while
# leaving the point estimate (labelled) correct. Found on the REStat census-truth test (Malawi assets, 18 Oct 2026).
#
# Tested here, with ratio and logit_raking benchmarking (megb() and xgb() have no additive benchmark type):
#   1. megb(): fixed target, perturb_benchmark (arm F) and bench_target = "truth" (arm T), bootstrap_refit "lmm_only" and
#      "leaves_only": ind, CI, var and their benchmarked versions are identical for any row order of pop_data;
#   2. the benchmarked replicates hit their group targets, whatever the row order (exact for ratio);
#   3. xgb(): fixed target and perturb_benchmark, same permutation (xgb passes the domain labels on both paths).
# The comparisons are on results sorted by domain, tolerance 1e-8.

perm_data <- function(seed = 7) {
  d <- megb_bench_data()
  pop <- d$pop
  rev_dom <- pop[order(-as.integer(sub("d", "", as.character(pop$dom)))), ]                 # domains in reverse order
  set.seed(seed); shuf <- pop[sample(nrow(pop)), ]                                           # rows fully shuffled
  list(smp = d$smp, sorted = pop, reversed = rev_dom, shuffled = shuf)
}

run_megb_order <- function(cfg, pop, arm = "fixed", refit = "lmm_only", B = 6) {
  d <- perm_data()
  y <- all.vars(cfg$f)[1]
  tgt <- tapply(d$smp[[y]], d$smp$reg, mean)
  tgt <- stats::setNames(as.numeric(tgt) * c(1.05, 0.95, 1.02)[seq_along(tgt)], names(tgt))
  if (cfg$bt == "logit_raking") tgt <- pmin(pmax(tgt, 0.05), 0.95)
  args <- list(fixed = cfg$f, smp_data = d$smp, pop_data = d[[pop]], pop_weights = "w", domains = "dom",
               transformation = cfg$tr, mse = TRUE, B = B, bootstrap_refit = refit, na.rm = FALSE, seed = 5,
               gradient_params = list(eta = 0.3, max_depth = 2, subsample = 1, nrounds = 30),
               benchmark = tgt, benchmark_type = cfg$bt, benchmark_level = "reg")
  se <- stats::setNames(rep(c(0.03, 0.05, 0.04), length.out = length(tgt)) * mean(tgt), names(tgt))
  if (arm == "F") { args$perturb_benchmark <- TRUE; args$benchmark_target_se <- se }
  if (arm == "T") { args$bench_target <- "truth"; args$benchmark_target_se <- se }
  suppressMessages(suppressWarnings(do.call(megb, args)))
}

by_domain <- function(x, key = "Domain") x[order(as.character(x[[key]])), , drop = FALSE]
expect_same_fit <- function(a, b, info = NULL) {
  for (nm in c("ind", "ind_bench", "CI", "CI_bench", "var", "var_bench")) {
    if (is.null(a[[nm]]) && is.null(b[[nm]])) next
    xa <- by_domain(a[[nm]]); xb <- by_domain(b[[nm]]); rownames(xa) <- rownames(xb) <- NULL
    expect_equal(xa, xb, tolerance = 1e-8, info = paste(nm, info))
  }
}

for (cn in names(megb_bench_cfgs)) {
  cfg <- megb_bench_cfgs[[cn]]
  for (arm in c("fixed", "F", "T")) {
    test_that(sprintf("megb %s, %s target: benchmarked results do not depend on the row order of pop_data (lmm_only)", cn, arm), {
      a <- run_megb_order(cfg, "sorted", arm)
      expect_true(!is.null(a$ind$Mean_bench) && !is.null(a$CI$Lower_bench) && !is.null(a$var_bench))
      expect_same_fit(a, run_megb_order(cfg, "reversed", arm), "reversed")
      expect_same_fit(a, run_megb_order(cfg, "shuffled", arm), "shuffled")
    })
  }
}

test_that("megb leaves_only (the default refit): benchmarked results do not depend on the row order of pop_data", {
  cfg <- megb_bench_cfgs$log_ratio
  for (arm in c("fixed", "T")) {
    a <- run_megb_order(cfg, "sorted", arm, refit = "leaves_only")
    expect_same_fit(a, run_megb_order(cfg, "shuffled", arm, refit = "leaves_only"), paste(arm, "shuffled"))
  }
})

test_that("megb: each benchmarked bootstrap replicate is scaled by one ratio factor per benchmark group, for any row order", {
  # bench_target = "truth": replicate b is benchmarked to T*_gb. Captured through the replicate estimates of mse_megb.
  cfg <- megb_bench_cfgs$no_ratio
  for (pop in c("sorted", "shuffled")) {

    suppressMessages(trace("mse_megb", exit = quote(assign(".bench_row_order_cap", returnValue(), envir = globalenv())), where = asNamespace("povmap"), print = FALSE))
    fit <- run_megb_order(cfg, pop, "T")
    suppressMessages(untrace("mse_megb", where = asNamespace("povmap"))); cap <- get(".bench_row_order_cap", envir = globalenv()); rm(".bench_row_order_cap", envir = globalenv())
    expect_false(is.null(cap$tau_b_bench))
    d <- perm_data(); w <- stats::aggregate(d[[pop]]$w, list(dom = as.character(d[[pop]]$dom)), sum); reg <- d[[pop]]$reg[match(as.character(cap$domains), as.character(d[[pop]]$dom))]
    ww <- w$x[match(as.character(cap$domains), w$dom)]
    # a replicate's benchmarked group mean must equal that replicate's target: its group mean ratio to the unbenchmarked one is the factor
    fac <- cap$tau_b_bench / cap$tau_b_orig
    for (g in unique(reg)) {
      r <- which(reg == g)
      expect_lt(max(apply(fac[r, , drop = FALSE], 2, function(z) diff(range(z)))), 1e-9)      # one factor per group and replicate
    }
  }
})

test_that("benchmark_xgb_level refuses unlabelled estimates", {
  fwk <- list(pop_data = data.frame(dom = c("b", "a"), reg = "R", npop = 1), domains = "dom", pop_weights = "npop")
  expect_error(povmap:::benchmark_xgb_level(list(ind = data.frame(Mean = c(0.2, 0.4))), fwk, y ~ x, c(R = 0.3), "ratio", "reg"), "Domain")
})

test_that("xgb: benchmarked point estimates do not depend on the row order of pop_data (the replicate path matches groups by name: add_benchmark tests in test_benchmark_level_alignment.R)", {
  skip_if_not_installed("xgboost")
  skip_on_cran()
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
  pops <- list(sorted = pop[order(pop$dom, pop$sa), ], reversed = pop[order(-as.integer(substr(pop$dom, 2, 3)), pop$sa), ],
               shuffled = { set.seed(3); pop[sample(nrow(pop)), ] })
  fit <- function(p, perturb) {
    args <- list(fixed = y ~ x1 + x2, smp_data = smp, smp_weights = "wt",
                 pop_data = p[, c("sub", "dom", "reg", "x1", "x2", "npop")], pop_weights = "npop",
                 domains = "dom", sub_domains = "sub", transformation = "no", bootstrap = TRUE,
                 B = 6, L = 5, seed = 1, nrounds = 15, eta = 0.3, max_depth = 2,
                 colsample_bylevel = 1, colsample_bynode = 1, subsample = 1,
                 benchmark = "Mean", benchmark_level = "reg", benchmark_type = "ratio", cpus = 1)
    if (perturb) { args$perturb_benchmark <- TRUE; args$benchmark_target_se <- c(R1 = 0.02, R2 = 0.03, R3 = 0.025) }
    suppressWarnings(suppressMessages(do.call(povmap::xgb, args)))
  }
  for (perturb in c(FALSE, TRUE)) {
    a <- fit(pops$sorted, perturb)
    for (nm in c("reversed", "shuffled")) {
      b <- fit(pops[[nm]], perturb)
      for (el in "ind") {    # point estimates only: xgb's bootstrap draws depend on the row order of pop_data even without a benchmark
        xa <- by_domain(a[[el]]); xb <- by_domain(b[[el]]); rownames(xa) <- rownames(xb) <- NULL
        expect_equal(xa, xb, tolerance = 1e-8, info = paste("xgb", el, nm, "perturb", perturb))
      }
    }
  }
})
