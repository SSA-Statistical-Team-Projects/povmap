# Shared by test-megb-perturb-benchmark.R and fixtures/make_megb_bench_fixture.R: a small megb problem with a
# benchmark level (three regions), numeric benchmark targets and bootstrap_refit = "lmm_only".
megb_bench_data <- function() {
  set.seed(7); D <- 14; n_pop <- 15; n_smp <- 6
  pop <- data.frame(dom = factor(rep(sprintf("d%02d", 1:D), each = n_pop)), x1 = rnorm(D * n_pop), x2 = runif(D * n_pop), w = rpois(D * n_pop, 20) + 1)
  pop$reg <- rep(c("A", "B", "C"), length.out = D)[as.integer(pop$dom)]
  u <- rnorm(D, 0, 0.4)[as.integer(pop$dom)]
  pop$y <- exp(0.4 * pop$x1 + u + rnorm(nrow(pop), 0, 0.4))
  pop$p <- plogis(0.5 * pop$x1 + u + rnorm(nrow(pop), 0, 0.3))
  pop$z <- 5 + 0.4 * pop$x1 + u + rnorm(nrow(pop), 0, 0.4)
  smp <- do.call(rbind, lapply(split(pop, pop$dom)[1:11], function(d) d[sample(nrow(d), n_smp), ]))
  list(smp = smp, pop = pop[, c("dom", "x1", "x2", "w", "reg")])
}

# cfg: list(f = formula, tr = transformation, bt = benchmark_type); `...` goes to megb() and overrides the defaults.
megb_bench_run <- function(cfg, ...) {
  d <- megb_bench_data()
  y <- all.vars(cfg$f)[1]
  tgt <- tapply(d$smp[[y]], d$smp$reg, mean)
  tgt <- stats::setNames(as.numeric(tgt) * c(1.05, 0.95, 1.02)[seq_along(tgt)], names(tgt))
  if (cfg$bt == "logit_raking") tgt <- pmin(pmax(tgt, 0.05), 0.95)
  args <- list(fixed = cfg$f, smp_data = d$smp, pop_data = d$pop, pop_weights = "w", domains = "dom",
               transformation = cfg$tr, mse = TRUE, B = 6, bootstrap_refit = "lmm_only", na.rm = FALSE, seed = 5,
               gradient_params = list(eta = 0.3, max_depth = 2, subsample = 1, nrounds = 30),
               benchmark = tgt, benchmark_type = cfg$bt, benchmark_level = "reg")
  extra <- list(...)
  args[names(extra)] <- extra
  suppressMessages(suppressWarnings(do.call(megb, args)))
}

megb_bench_cfgs <- list(
  arcsin_logit = list(f = p ~ x1 + x2, tr = "arcsin", bt = "logit_raking"),
  log_ratio    = list(f = y ~ x1 + x2, tr = "log", bt = "ratio"),
  no_ratio     = list(f = z ~ x1 + x2, tr = "no", bt = "ratio"))
