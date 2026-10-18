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

# Run megb() and capture what mse_megb() returns (the B x areas replicate matrices) and the calls of the benchmark function.
capture_megb <- function(cfg, pop, arm, refit) {
  ns <- asNamespace("povmap"); assign(".bro_calls", list(), envir = globalenv())
  suppressMessages({
    trace("mse_megb", exit = quote(assign(".bro_cap", returnValue(), envir = globalenv())), where = ns, print = FALSE)
    trace("benchmark_xgb_level", quote(assign(".bro_calls", c(get(".bro_calls", envir = globalenv()),
          list(data.frame(Domain = as.character(point_estim$ind$Domain), Mean = point_estim$ind$Mean))), envir = globalenv())), where = ns, print = FALSE)
  })
  on.exit(suppressMessages({ untrace("mse_megb", where = ns); untrace("benchmark_xgb_level", where = ns) }), add = TRUE)
  fit <- run_megb_order(cfg, pop, arm, refit = refit)
  cap <- get(".bro_cap", envir = globalenv()); calls <- get(".bro_calls", envir = globalenv())
  rm(list = c(".bro_cap", ".bro_calls"), envir = globalenv())
  list(fit = fit, cap = cap, calls = calls)
}

test_that("megb bench_target = truth: each replicate hits its target (the group truth plus that replicate's draw, clamped), for any row order", {
  # The target of replicate b in group g is the population-weighted group mean of tau_star_orig[, b] plus draw[b, g], clamped to
  # the bounds of the transformation. The weighted group mean of the benchmarked replicate tau_b_bench[, b] must equal it.
  for (cn in c("no_ratio", "arcsin_logit")) {
    cfg <- megb_bench_cfgs[[cn]]; bnd <- povmap:::.megb_bench_bounds(cfg$tr)
    for (pop in c("sorted", "shuffled")) for (refit in c("lmm_only", "leaves_only")) {
      r <- capture_megb(cfg, pop, "T", refit); cap <- r$cap; draws <- r$fit$benchmark_target_draws
      expect_false(is.null(cap$tau_b_bench)); expect_false(is.null(draws))
      d <- perm_data()[[pop]]; wt <- tapply(d$w, as.character(d$dom), sum); dm <- as.character(cap$domains)
      w <- as.numeric(wt[dm]); reg <- d$reg[match(dm, as.character(d$dom))]
      keep <- cap$bench_keep; if (is.null(keep)) keep <- seq_len(ncol(cap$tau_b_bench))
      expect_equal(length(keep), ncol(cap$tau_b_bench))
      wm <- function(M, rows, j) sum(M[rows, j] * w[rows]) / sum(w[rows])
      lg <- function(p) log(p / (1 - p))
      for (g in colnames(draws)) {
        rows <- which(reg == g)
        for (j in seq_along(keep)) {
          tgt <- min(max(wm(cap$tau_star_orig, rows, keep[j]) + draws[keep[j], g], bnd[["lo"]]), bnd[["hi"]])
          expect_equal(wm(cap$tau_b_bench, rows, j), tgt, tolerance = 1e-6, info = paste(cn, pop, refit, g, "replicate", keep[j]))
          if (cfg$bt == "logit_raking" && tgt > 1e-6 && tgt < 1 - 1e-6) {   # logit raking: one constant logit shift per group and replicate
            sh <- lg(cap$tau_b_bench[rows, j]) - lg(cap$tau_b_orig[rows, keep[j]])
            expect_lt(diff(range(sh)), 1e-5, label = paste(cn, pop, refit, g, "replicate", keep[j], "logit shift range"))
          }
        }
      }
    }
  }
})

test_that("megb leaves_only: the domain labels handed to the benchmark are in the order of the replicate's estimates (mean_preds$dom_name vs mean_boot_orig)", {
  # Each replicate's benchmark_xgb_level() call carries (Domain, Mean): the Mean of a domain must be that domain's unbenchmarked
  # replicate estimate in tau_b_orig, so the labels are in the order of mean_boot_orig whatever the row order of pop_data.
  cfg <- megb_bench_cfgs$no_ratio
  for (pop in c("sorted", "reversed", "shuffled")) {
    r <- capture_megb(cfg, pop, "T", "leaves_only"); cap <- r$cap; dm <- as.character(cap$domains)
    rc <- r$calls[-1]                                           # the first call is the point estimate
    expect_gte(length(rc), 1)
    for (cl in rc) {
      ok <- vapply(seq_len(ncol(cap$tau_b_orig)), function(j) isTRUE(all.equal(cl$Mean, cap$tau_b_orig[match(cl$Domain, dm), j], tolerance = 1e-10)), NA)
      expect_true(any(ok), info = paste(pop, "a replicate's labelled estimates match no column of tau_b_orig"))
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
