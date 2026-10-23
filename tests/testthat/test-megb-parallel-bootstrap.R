# bootstrap_cores > 1 with bootstrap_refit = "leaves_only" runs the replicates on a PSOCK cluster. The replicate draws no random numbers, so the
# result must be identical (not just close) to the sequential loop for the same seed. Needs the installed package on the workers.
test_that("parallel leaves_only bootstrap is identical to the sequential loop", {
  skip_on_cran(); skip_if_not_installed("xgboost"); skip_if_not_installed("collapse"); skip_if_not_installed("jsonlite")
  skip_if_not("povmap" %in% rownames(utils::installed.packages()))
  set.seed(42)
  D <- 40; n_pop <- 30; n_smp <- 6
  pop <- data.frame(dom = factor(rep(sprintf("d%02d", 1:D), each = n_pop)),
                    x1 = rnorm(D * n_pop), x2 = runif(D * n_pop), w = rpois(D * n_pop, 20) + 1)
  u <- rnorm(D, 0, 0.3)
  pop$y <- 1 + 0.8 * pop$x1 - 0.5 * pop$x2^2 + u[as.integer(pop$dom)] + rnorm(nrow(pop), 0, 0.5)
  smp <- do.call(rbind, lapply(split(pop, pop$dom)[1:30], function(d) d[sample(nrow(d), n_smp), ]))
  smp$wt <- runif(nrow(smp), 0.5, 2)
  gp <- list(eta = 0.1, max_depth = 3, subsample = 0.8, colsample_bytree = 0.9, nrounds = 60, nthread = 1)
  old <- Sys.getenv("OMP_NUM_THREADS"); Sys.setenv(OMP_NUM_THREADS = "1"); on.exit(Sys.setenv(OMP_NUM_THREADS = old), add = TRUE)
  fit <- function(cores, leaf, tr) {
    o <- options(povmap.leaves_only.leaf_aggregate = leaf); on.exit(options(o))
    suppressMessages(suppressWarnings(povmap::megb(
      fixed = y ~ x1 + x2, smp_data = smp, smp_weights = "wt", pop_data = pop[, c("dom", "x1", "x2", "w")],
      pop_weights = "w", domains = "dom", transformation = tr, mse = TRUE, B = 6, bootstrap_refit = "leaves_only",
      bootstrap_cores = cores, na.rm = FALSE, seed = 11, gradient_params = gp)))
  }
  for (cfg in list(list(TRUE, "no"), list(FALSE, "no"))) {   # leaf aggregation and per-cell prediction
    a <- fit(0, cfg[[1]], cfg[[2]]); b <- fit(2, cfg[[1]], cfg[[2]])
    for (nm in c("ind", "MSE", "var", "CI")) expect_identical(a[[nm]], b[[nm]], info = paste(nm, cfg[[1]], cfg[[2]]))
  }
})

# The configuration of the Nigeria / DRC runs: arcsin with the plain back-transform (no leaf shortcut: per-cell predictions), logit raking to group targets,
# bench_target = "truth", boot_area_var = "centred" and boot_synthetic = TRUE.
test_that("parallel leaves_only bootstrap is identical for arcsin plain, logit raking, truth target, centred area variance and boot_synthetic", {
  skip_on_cran(); skip_if_not_installed("xgboost"); skip_if_not_installed("collapse"); skip_if_not_installed("jsonlite")
  skip_if_not("povmap" %in% rownames(utils::installed.packages()))
  set.seed(43)
  D <- 40; n_pop <- 30; n_smp <- 6
  dom <- factor(rep(sprintf("d%02d", 1:D), each = n_pop))
  pop <- data.frame(dom = dom, x1 = rnorm(D * n_pop), x2 = runif(D * n_pop), w = rpois(D * n_pop, 20) + 1)
  pop$grp <- paste0("g", (as.integer(pop$dom) - 1) %/% 10 + 1)                       # four benchmark groups of ten domains
  u <- rnorm(D, 0, 0.3)
  pop$y <- stats::plogis(-0.3 + 0.7 * pop$x1 - 0.5 * pop$x2 + u[as.integer(pop$dom)] + rnorm(nrow(pop), 0, 0.3))
  sampled <- sprintf("d%02d", c(1:8, 11:18, 21:28, 31:38))                           # eight unsampled domains, two in each group
  smp <- do.call(rbind, lapply(split(pop, pop$dom)[sampled], function(d) d[sample(nrow(d), n_smp), ]))
  smp$wt <- runif(nrow(smp), 0.5, 2)
  bm <- vapply(split(smp, smp$grp), function(d) stats::weighted.mean(d$y, d$wt), 0)
  gp <- list(eta = 0.1, max_depth = 3, subsample = 0.8, colsample_bytree = 0.9, nrounds = 60, nthread = 1)
  old <- Sys.getenv("OMP_NUM_THREADS"); Sys.setenv(OMP_NUM_THREADS = "1"); on.exit(Sys.setenv(OMP_NUM_THREADS = old), add = TRUE)
  fit <- function(cores) suppressMessages(suppressWarnings(povmap::megb(
    fixed = y ~ x1 + x2, smp_data = smp, smp_weights = "wt", pop_data = pop[, c("dom", "x1", "x2", "w", "grp")], pop_weights = "w", domains = "dom",
    transformation = "arcsin", arcsin_backtransform = "plain", mse = TRUE, B = 6, bootstrap_refit = "leaves_only", bootstrap_cores = cores, na.rm = FALSE,
    seed = 11, gradient_params = gp, benchmark = bm, benchmark_type = "logit_raking", benchmark_level = "grp", bench_target = "truth", perturb_benchmark = FALSE,
    benchmark_target_se = setNames(rep(0.03, length(bm)), names(bm)), boot_area_var = "centred", boot_synthetic = TRUE)))
  a <- fit(0); b <- fit(2)
  for (nm in c("ind", "ind_bench", "var", "var_bench", "CI", "CI_bench", "boot_diag")) expect_identical(a[[nm]], b[[nm]], info = nm)
  expect_false(is.null(a$boot_diag$synthetic))                                        # boot_synthetic was in effect
})

test_that("the parallel bootstrap reports failed replicates and refuses a worker with a different build", {
  skip_on_cran(); skip_if_not("povmap" %in% rownames(utils::installed.packages()))
  # a worker that loaded another povmap / xgboost than the main process stops the call (checked without a cluster: the worker-side check is a plain function)
  id <- povmap:::.megb_pkg_id()
  expect_null(povmap:::.megb_par_init(NULL, id))
  expect_error(povmap:::.megb_par_init(NULL, modifyList(id, list(povmap_version = "0.0.0"))), "different build")
  expect_error(povmap:::.megb_par_init(NULL, modifyList(id, list(remote_sha = "deadbeef"))), "different build")
  expect_error(povmap:::.megb_par_init(NULL, modifyList(id, list(povmap_lib = "Z:/elsewhere"))), "different build")
  # replicates that come back empty are counted and their messages returned in a warning (here every replicate fails: the state sent to the worker is empty)
  old <- Sys.getenv("OMP_NUM_THREADS"); Sys.setenv(OMP_NUM_THREADS = "1"); on.exit(Sys.setenv(OMP_NUM_THREADS = old), add = TRUE)
  S <- list(y_star_smp = matrix(1, 3, 3), opt_em = FALSE, MaxIterations = 1L)
  expect_warning(res <- povmap:::.megb_leaf_rep_parallel(3L, S, 1L), "3 of 3 replicates failed on the workers.*leaf EM failed at iter")
  expect_length(res, 3L)
  expect_true(all(vapply(res, function(r) is.null(r$Mean_boot) && is.null(r$fail), NA)))   # the failure text is not left in the replicate lists
})
