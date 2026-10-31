# Data and a megb() call for test-megb-boot-checkpoint.R (also sourced by the child R process of the kill test): the configuration of the Nigeria / DRC
# runs (arcsin plain, logit raking to group targets, truth target, centred area variance, boot_synthetic), so that the fingerprint covers a benchmark closure.
ck_data <- function() {
  set.seed(43)
  D <- 40; n_pop <- 30; n_smp <- 6
  dom <- factor(rep(sprintf("d%02d", 1:D), each = n_pop))
  pop <- data.frame(dom = dom, x1 = rnorm(D * n_pop), x2 = runif(D * n_pop), w = rpois(D * n_pop, 20) + 1)
  pop$grp <- paste0("g", (as.integer(pop$dom) - 1) %/% 10 + 1)
  u <- rnorm(D, 0, 0.3)
  pop$y <- stats::plogis(-0.3 + 0.7 * pop$x1 - 0.5 * pop$x2 + u[as.integer(pop$dom)] + rnorm(nrow(pop), 0, 0.3))
  sampled <- sprintf("d%02d", c(1:8, 11:18, 21:28, 31:38))
  smp <- do.call(rbind, lapply(split(pop, pop$dom)[sampled], function(d) d[sample(nrow(d), n_smp), ]))
  smp$wt <- runif(nrow(smp), 0.5, 2)
  bm <- vapply(split(smp, smp$grp), function(d) stats::weighted.mean(d$y, d$wt), 0)
  list(smp = smp, pop = pop, bm = bm)
}

ck_fit <- function(cores = 0, B = 50, d = ck_data(), seed = 11, ...) {
  gp <- list(eta = 0.1, max_depth = 3, subsample = 0.8, colsample_bytree = 0.9, nrounds = 60, nthread = 1)
  suppressWarnings(povmap::megb(
    fixed = y ~ x1 + x2, smp_data = d$smp, smp_weights = "wt", pop_data = d$pop[, c("dom", "x1", "x2", "w", "grp")], pop_weights = "w", domains = "dom",
    transformation = "arcsin", arcsin_backtransform = "plain", mse = TRUE, B = B, bootstrap_refit = "leaves_only", bootstrap_cores = cores, na.rm = FALSE,
    seed = seed, gradient_params = gp, benchmark = d$bm, benchmark_type = "logit_raking", benchmark_level = "grp", bench_target = "truth", perturb_benchmark = FALSE,
    benchmark_target_se = setNames(rep(0.03, length(d$bm)), names(d$bm)), boot_area_var = "centred", boot_synthetic = TRUE, ...))
}

ck_same <- function(a, b, info = "") {
  for (nm in c("ind", "ind_bench", "MSE", "var", "var_bench", "CI", "CI_bench", "boot_diag"))
    testthat::expect_identical(a[[nm]], b[[nm]], info = paste(nm, info))
}
