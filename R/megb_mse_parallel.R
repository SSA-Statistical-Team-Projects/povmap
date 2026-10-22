# Parallel leaves_only bootstrap (DRAFT, branch megb-speed-bootstrap off 2913cc9; not on any shared branch or tag).
#
# The leaves_only replicate loop of mse_megb() draws no random numbers (y_star_smp, the bootstrap truth and the benchmark target draws
# are all made before the loop), and every replicate depends only on its column of y_star_smp and on objects built once per fit. The
# replicates are therefore independent and give the same numbers in any order on any worker.
#
# .megb_leaf_rep() is a copy of the loop body in mse_megb() (the sequential loop there is untouched and stays the default), with
# `next` turned into return() and the options read in the main process passed in S. mse_megb() calls .megb_leaf_rep_parallel()
# when bootstrap_refit = "leaves_only", bootstrap_cores > 1 and no diagnostic environment is set. Bit-identity of the two paths is
# tested in tests/testthat/test-megb-parallel-bootstrap.R and by the script that comes with the draft.

.megb_par_env <- new.env(parent = emptyenv())

# domain mean of a cell-level vector (the mse_megb() agg_dom closure, standalone so that no function frame is serialised)
.megb_agg_dom <- function(v, pop_dom_idx, pw_pop) {
  if (is.null(pw_pop)) as.numeric(tapply(v, pop_dom_idx, mean))
  else as.numeric(tapply(v * pw_pop, pop_dom_idx, sum) / tapply(pw_pop, pop_dom_idx, sum))
}

.megb_leaf_rep <- function(i, S) {
  CFB <- S$CFB; X_smp_mat <- S$X_smp_mat; X_pop_mat <- S$X_pop_mat; y_star_smp <- S$y_star_smp
  w_ref <- S$w_ref; fit_data <- S$fit_data; lmm_formula <- S$lmm_formula; smp_weights_vec <- S$smp_weights_vec
  .use_leaf <- S$use_leaf; .cfp <- S$cfp; .cfu <- S$cfu; MaxIterations <- S$MaxIterations; ErrorTolerance <- S$ErrorTolerance
  unit_pred_smp <- S$unit_pred_smp; gb_smp <- S$gb_smp; orig_booster <- S$orig_booster; orig_nrounds <- S$orig_nrounds
  LA <- S$LA; newdat_d <- S$newdat_d; pop_dom_idx <- S$pop_dom_idx; pw_pop <- S$pw_pop; corrected_bt <- S$corrected_bt
  domains <- S$domains; boot_synthetic <- S$boot_synthetic; benchmark_fn <- S$benchmark_fn; bm_target <- S$bm_target; smp_data <- S$smp_data
  agg_dom <- function(v) .megb_agg_dom(v, pop_dom_idx, pw_pop)

  y_i <- y_star_smp[, i]
  em  <- tryCatch(
    .megb_leaf_em(CFB, X_smp_mat, X_pop_mat, y_i, w_ref, fit_data, lmm_formula, smp_weights_vec,
                  use_leaf = .use_leaf, predict = .cfp, predict_unsampled = .cfu,
                  max_iterations = if (isFALSE(S$opt_em)) 0L else MaxIterations,
                  error_tolerance = ErrorTolerance,
                  start_re = if (identical(S$opt_em_start, "point")) unit_pred_smp - gb_smp else NULL),
    error = function(e) { message("  leaf EM failed at iter ", i, ": ", conditionMessage(e)); NULL })
  if (is.null(em))
    return(list(Mean_boot = NULL, Mean_boot_orig = NULL, Mean_boot_bench = NULL,
                error_sd_boot = NA_real_, ran_eff_sd_boot = NA_real_, em_iterations = NA_integer_))
  gb_refreshed <- tryCatch(
    .megb_refresh(orig_booster, X_smp_mat, em$target, w_ref, orig_nrounds),
    error = function(e) { message("  refresh failed at iter ", i, ": ", conditionMessage(e)); NULL })
  if (is.null(gb_refreshed))
    return(list(Mean_boot = NULL, Mean_boot_orig = NULL, Mean_boot_bench = NULL,
                error_sd_boot = NA_real_, ran_eff_sd_boot = NA_real_))
  if (.use_leaf) {
    LM <- .megb_leaf_domain_means(LA, gb_refreshed)
    gb_pop_boot <- NULL
  } else {
    gb_pop_boot <- as.numeric(predict(gb_refreshed, X_pop_mat))
  }
  {
    cfr <- em$cfr
    if (.cfp && .use_leaf) {
      LM$unweighted[CFB$pos_all] <- cfr$unweighted
      LM$weighted[CFB$pos_all]   <- cfr$weighted
    } else if (.cfp) {
      gb_pop_boot[CFB$ip_all] <- cfr$gb_pop
    }
    if (.cfu && length(CFB$posU) && .use_leaf) {
      LM$unweighted[CFB$posU] <- cfr$u_unweighted
      LM$weighted[CFB$posU]   <- cfr$u_weighted
    } else if (.cfu && length(CFB$posU)) {
      gb_pop_boot[CFB$ipU] <- cfr$u_gb_pop
    }
  }
  lmer_boot <- em$lmer
  if (is.null(lmer_boot))
    return(list(Mean_boot = NULL, Mean_boot_orig = NULL, Mean_boot_bench = NULL,
                error_sd_boot = NA_real_, ran_eff_sd_boot = NA_real_))
  re_per_d <- stats::predict(lmer_boot, newdata = newdat_d, allow.new.levels = TRUE)
  fe_boot  <- as.numeric(lme4::fixef(lmer_boot))[1]
  u_d_boot <- as.numeric(re_per_d) - fe_boot
  gb_pop_d_mean_boot <- if (.use_leaf) LM$unweighted
                        else as.numeric(tapply(gb_pop_boot, pop_dom_idx, mean))
  mean_boot_t        <- gb_pop_d_mean_boot + u_d_boot
  synth_boot <- NULL
  if (isTRUE(boot_synthetic)) {
    synth_boot <- if (!is.null(corrected_bt) && .use_leaf) LM$weighted
                  else if (!is.null(corrected_bt))
                    agg_dom(corrected_bt(gb_pop_boot, .megb_re_var(corrected_bt, lmer_boot, rep(".synthetic.", length(domains)))[pop_dom_idx]))
                  else gb_pop_d_mean_boot
  }
  mean_boot_orig     <- if (!is.null(corrected_bt) && .use_leaf) {
    LM$weighted + u_d_boot
  } else if (!is.null(corrected_bt)) {
    cell_orig_boot <- corrected_bt(gb_pop_boot + u_d_boot[pop_dom_idx],
                                   .megb_re_var(corrected_bt, lmer_boot, domains)[pop_dom_idx])
    agg_dom(cell_orig_boot)
  } else mean_boot_t
  mean_boot_bench <- if (!is.null(benchmark_fn)) {
    bs_for_bench <- smp_data
    bs_for_bench$y_star <- y_star_smp[, i]
    if (!is.null(bm_target)) attr(bs_for_bench, "bm_target") <- bm_target[, i]
    benchmark_fn(mean_boot_orig, bs_for_bench, domains)
  } else NULL
  ran_sd <- as.data.frame(lme4::VarCorr(lmer_boot))$sdcor[1]
  err_sd <- stats::sigma(lmer_boot)
  out <- list(Mean_boot = mean_boot_t, Mean_boot_orig = mean_boot_orig, Mean_boot_bench = mean_boot_bench,
              error_sd_boot = err_sd, ran_eff_sd_boot = ran_sd, em_iterations = em$iterations)
  if (isTRUE(boot_synthetic)) out$Mean_boot_synth <- synth_boot
  out
}

.megb_par_init <- function(S) { assign("S", S, envir = .megb_par_env); invisible(NULL) }
.megb_par_run  <- function(i) .megb_leaf_rep(i, get("S", envir = .megb_par_env))

# Start the PSOCK cluster, retrying after a short pause: makeCluster() can fail on a busy machine when the port it picked is taken (seen once in testing).
.megb_make_cluster <- function(cores, tries = 3L, pause = 5) {
  err <- NULL
  for (a in seq_len(tries)) {
    cl <- tryCatch(parallel::makeCluster(cores), error = function(e) e)
    if (!inherits(cl, "error")) return(cl)
    err <- cl
    if (a < tries) { message("leaves_only bootstrap: cluster start failed (", conditionMessage(err), "); retrying in ", pause, " s"); Sys.sleep(pause) }
  }
  stop(err)
}

# Run the B replicates on `cores` PSOCK workers (each gets S once) and return the list in replicate order.
# The workers inherit the environment of this process (R_LIBS, R_PROFILE_USER, OMP_NUM_THREADS), so they load the same povmap and
# xgboost builds; set OMP_NUM_THREADS = 1 first for single-threaded refreshes, as the sequential default does.
.megb_leaf_rep_parallel <- function(B, S, cores) {
  if (isTRUE(getOption("povmap.leaves_only.par_diag", FALSE)))   # serialised size of each object every worker receives
    message("leaves_only bootstrap, MB sent to each worker: ", paste(names(S), round(vapply(S, function(x) length(serialize(x, NULL)), 0) / 1e6, 1), collapse = ", "))
  cl <- .megb_make_cluster(cores)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  parallel::clusterCall(cl, .megb_par_init, S)
  parallel::parLapplyLB(cl, seq_len(B), .megb_par_run, chunk.size = 2L)
}
