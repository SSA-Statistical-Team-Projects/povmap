# Internal: parametric bootstrap MSE for megb_em
#' @importFrom purrr map
#' @importFrom furrr future_map furrr_options
#' @importFrom future plan multisession multicore sequential
#' @importFrom dplyr group_by summarise

mse_megb <- function(Y, X, dom_name, smp_data, model, error_sd, pop_data,
                     B, initial_random_effects, ErrorTolerance, MaxIterations,
                     cov_names, gradient_params = list(), formula_random_effects,
                     bootstrap_cores, seed, gbm_engine,
                     unit_pred_smp, gb_smp = NULL, unit_preds,
                     bootstrap_refit = "leaves_only",
                     corrected_bt = NULL,
                     benchmark_fn = NULL,
                     smp_weights_col = NULL,
                     weightedBS = TRUE,
                     ...) {
  # smp_weights_col: name of a column on smp_data containing observation
  #   weights. Pulling weights from a column (rather than passing a vector)
  #   guarantees row alignment under any sorting / subsetting smp_data
  #   undergoes within this function. When NULL (or the column doesn't exist),
  #   the bootstrap operates as if all observations are equally weighted.
  # weightedBS: when TRUE and smp_weights_col is supplied, residual sampling
  #   (both level-1 in block_sample and level-2 for u_d* draws) uses
  #   probabilities proportional to weights. Mirrors xgb's weightedBS arg.

  # smp_weights_vec is extracted from the column AFTER the domain sort below
  # to guarantee alignment with smp_data row order.

  # Function-scoped definition so any branch that references collect_diag has
  # it in scope; branches that don't use it just leave it FALSE.
  collect_diag <- !is.null(getOption("povmap.leaves_only.diag_env"))

  # The bootstrap is an EBLUP-style parametric residual bootstrap, refit on
  # the *original* smp_data (with bootstrapped y) — not on a pop_data subsample.
  # Refitting on a pop subsample (the previous behaviour) made each bootstrap
  # GB train on a different X distribution than the original GB and pushed
  # error_sd_boot far below the original error_sd, inflating prediction MSE.

  # Sort smp_data and pop_data by domain so the per-domain block samplers,
  # tapply aggregation, and sample_ui all operate on consistent orderings.
  sort_smp <- order(as.character(smp_data[[dom_name]]))
  smp_data <- smp_data[sort_smp, , drop = FALSE]
  Y        <- Y[sort_smp]
  if (!is.null(unit_pred_smp))   unit_pred_smp   <- unit_pred_smp[sort_smp]
  if (!is.null(gb_smp))          gb_smp          <- gb_smp[sort_smp]
  # Re-extract weights from the (now-sorted) column to guarantee alignment.
  smp_weights_vec <- if (!is.null(smp_weights_col) &&
                         smp_weights_col %in% colnames(smp_data)) {
    as.numeric(smp_data[[smp_weights_col]])
  } else NULL

  # Defensive: if for any reason the extracted vector still doesn't match the
  # smp_data row count, disable weighting rather than crash.
  if (!is.null(smp_weights_vec) && length(smp_weights_vec) != nrow(smp_data)) {
    warning("smp_weights_col '", smp_weights_col, "' length (",
            length(smp_weights_vec), ") != nrow(smp_data) (",
            nrow(smp_data), "). Disabling weighted bootstrap.")
    smp_weights_vec <- NULL
  }

  # Rescale weights WITHIN each domain (Pfeffermann scaling), matching the
  # original fit in em_gb_lmm so the bootstrap LMM refits and weighted residual
  # sampling use the same weight treatment. Weights sum to the domain sample
  # size (average 1) within each domain. smp_data is already domain-sorted above.
  if (!is.null(smp_weights_vec)) {
    if (any(is.na(smp_weights_vec))) {
      warning("smp_weights had ", sum(is.na(smp_weights_vec)),
              " NA values; replacing with 0 before rescaling.")
      smp_weights_vec[is.na(smp_weights_vec)] <- 0
    }
    dom  <- as.character(smp_data[[dom_name]])
    dsum <- ave(smp_weights_vec, dom, FUN = function(z) sum(z, na.rm = TRUE))
    dn   <- ave(smp_weights_vec, dom, FUN = length)
    keep <- is.finite(dsum) & dsum > 0
    smp_weights_vec[keep]  <- smp_weights_vec[keep] * dn[keep] / dsum[keep]
    smp_weights_vec[!keep] <- 1
  }

  sort_pop   <- order(as.character(pop_data[[dom_name]]))
  pop_data   <- pop_data[sort_pop, , drop = FALSE]
  unit_preds <- unit_preds[sort_pop, , drop = FALSE]

  # Fixed-effect-only predictors (GB output, *not* GB + BLUP). Falling back to
  # the combined predictions would re-introduce the area-effect double counting
  # the comment chain in the prior version warned about.
  if (!is.null(unit_preds$gb_pop)) {
    gb_pop_vec <- unit_preds$gb_pop
  } else {
    warning("unit_preds lacks a 'gb_pop' column; falling back to combined ",
            "unit predictions. Bootstrap MSE will be inflated.")
    gb_pop_vec <- unit_preds$unit_preds
  }
  if (is.null(gb_smp)) {
    warning("gb_smp not supplied; falling back to combined unit_pred_smp. ",
            "Bootstrap will inject a doubled area effect into y_star_smp.")
    gb_smp <- unit_pred_smp
  }

  domains  <- as.character(unique(pop_data[[dom_name]]))     # alphabetical after sort
  in_samp  <- domains %in% as.character(unique(smp_data[[dom_name]]))
  n_d_pop  <- as.numeric(table(pop_data[[dom_name]]))        # row counts in domain order
  n_d_smp_tab <- table(smp_data[[dom_name]])
  # n_d_smp aligned to `domains`: zero for OOS domains
  n_d_smp <- as.numeric(n_d_smp_tab[match(domains, names(n_d_smp_tab))])
  n_d_smp[is.na(n_d_smp)] <- 0

  ran_obj  <- ran_comp(
    Y = Y, smp_data = smp_data, unit_pred_smp = unit_pred_smp,
    error_sd = error_sd, dom_name = dom_name,
    cov_names = cov_names, model = model,
    weights = smp_weights_vec
  )
  ran_effs <- ran_obj$ran_effs
  gb_res   <- ran_obj$gb_res
  smp_data <- ran_obj$smp_data
  area_w   <- ran_obj$area_w  # per-domain summed weights (or NULL)

  # ── y_star at SMP level: y_star_smp = gb_smp + u_d*[d, b] + e_ij*[i, b] ──
  # u_d* is sampled with replacement once per area per iteration, then
  # broadcast to each smp unit in that domain. e_ij* is block-sampled from
  # gb_res at smp size per domain.
  smp_pred_mat <- matrix(gb_smp, nrow = length(gb_smp), ncol = B)

  # Block sample residuals at SMP size: pass smp_data as the size reference
  # so block_sample yields sum(smp_data[d] == d) draws per in-sample domain.
  block_sample_e_smp <- function(x) {
    block_sample(domains = domains, in_samp = in_samp, smp_data = smp_data,
                 dom_name = dom_name, pop_data = smp_data, gb_res = gb_res,
                 weights = smp_weights_vec, weightedBS = weightedBS)
  }

  # u_d* sampler: one fresh draw per domain per iteration, length D.
  # When weightedBS + smp_weights_vec are supplied, sample with probability
  # proportional to area-summed weights (matching xgb's pop_area_d treatment).
  ud_prob <- if (!is.null(area_w) && isTRUE(weightedBS)) area_w / sum(area_w) else NULL
  sample_ud <- function(x) sample(ran_effs, size = length(domains), replace = TRUE,
                                  prob = ud_prob)

  e_smp    <- apply(matrix(NA, nrow = length(gb_smp), ncol = B), 2, block_sample_e_smp)
  u_d_star <- apply(matrix(NA, nrow = length(domains), ncol = B), 2, sample_ud)

  # Map u_d* to smp-row level via domain index.
  smp_dom_idx <- match(as.character(smp_data[[dom_name]]), domains)
  u_smp       <- u_d_star[smp_dom_idx, , drop = FALSE]
  smp_data$gb_res <- NULL

  y_star_smp <- smp_pred_mat + u_smp + e_smp

  # ── tau_star at POP level (analytical) ──
  # tau_star_d[b] = mean over pop_d of (gb_pop + u_d*[d, b]) ≈ gb_pop_d_mean + u_d*[d, b]
  # The within-domain mean of e_ij* averages to ≈ 0 for typical pop sizes, so
  # we don't materialise pop-level e — saves memory and avoids extra noise.
  pop_dom_idx   <- match(as.character(pop_data[[dom_name]]), domains)
  gb_pop_d_mean <- as.numeric(tapply(gb_pop_vec, pop_dom_idx, mean))
  tau_star      <- gb_pop_d_mean + u_d_star    # broadcasts: D × B

  # ── Bootstrap "samples" are the original smp_data with y_star injected ──
  # No subsampling: each iter sees the same smp X distribution as the original
  # fit, which keeps gb_pop_boot ≈ gb_pop on average and (crucially) keeps
  # error_sd_boot in line with the original error_sd.
  boots_sample <- vector(mode = "list", length = B)
  for (i in seq_len(B)) {
    bs                <- smp_data
    bs$y_star         <- y_star_smp[, i]
    boots_sample[[i]] <- bs
  }

  # ── Branch: lmm_only / leaves_only / full bootstrap refit ───────────────────
  # bootstrap_refit = "lmm_only": treat the gradient booster as fixed
  # across iterations and refit only the linear mixed model on bootstrap
  # residuals r = u_d* + e_smp* = y_star_smp - gb_smp. This is the standard
  # EBLUP parametric bootstrap (Prasad–Rao, Hall–Maiti, Pfeffermann–Tiller),
  # which yields the textbook leading-order MSE g_1 = γ_d σ_e² / n_d. It
  # isolates random-effect uncertainty, runs an order of magnitude faster than
  # refitting the booster each iteration, and avoids contaminating MSE with
  # GB-fit drift.
  if (bootstrap_refit == "lmm_only") {
    lmm_formula <- stats::as.formula(paste0("r ~ 1 + ", formula_random_effects))
    # newdata for predicting per-domain BLUPs (one row per domain, in `domains`
    # order). predict() with allow.new.levels=TRUE returns the fixed-effect
    # intercept for OOS levels, so re_per_d - fe = u_d_boot (zero for OOS).
    newdat_d <- stats::setNames(
      data.frame(domains, stringsAsFactors = FALSE), dom_name
    )
    if (is.factor(smp_data[[dom_name]])) {
      newdat_d[[dom_name]] <- factor(domains, levels = levels(smp_data[[dom_name]]))
    }
    boots_models <- vector("list", B)
    fit_data     <- smp_data       # we just append/overwrite the `r` column
    for (i in seq_len(B)) {
      if (i %% max(1L, B %/% 10L) == 0L)
        message("LMM-only bootstrap iteration ", i, " of ", B)

      # Residual = y_star_smp - gb_smp = u_smp[, i] + e_smp[, i] (saves an add).
      fit_data$r <- u_smp[, i] + e_smp[, i]
      if (!is.null(smp_weights_vec)) fit_data$.megb_w <- smp_weights_vec

      lmer_boot <- tryCatch(
        suppressMessages(suppressWarnings(
          lme4::lmer(
            lmm_formula, data = fit_data, REML = TRUE,
            weights = if (!is.null(smp_weights_vec)) fit_data$.megb_w else NULL
          )
        )),
        error = function(e) NULL
      )
      if (is.null(lmer_boot)) {
        boots_models[[i]] <- list(
          Mean_boot = NULL, Mean_boot_orig = NULL, Mean_boot_bench = NULL,
          error_sd_boot = NA_real_, ran_eff_sd_boot = NA_real_
        )
        next
      }

      re_per_d <- stats::predict(lmer_boot, newdata = newdat_d,
                                 allow.new.levels = TRUE)
      fe_boot  <- as.numeric(lme4::fixef(lmer_boot))[1]   # intercept
      u_d_boot <- as.numeric(re_per_d) - fe_boot   # length D, 0 for OOS

      mean_boot_t    <- gb_pop_d_mean + u_d_boot
      mean_boot_orig <- if (!is.null(corrected_bt)) corrected_bt(mean_boot_t)
                        else mean_boot_t

      mean_boot_bench <- if (!is.null(benchmark_fn)) {
        bs_for_bench <- smp_data
        bs_for_bench$y_star <- gb_smp + fit_data$r   # = y_star_smp[, i]
        benchmark_fn(mean_boot_orig, bs_for_bench)
      } else NULL

      ran_sd <- as.data.frame(lme4::VarCorr(lmer_boot))$sdcor[1]
      err_sd <- stats::sigma(lmer_boot)

      boots_models[[i]] <- list(
        Mean_boot       = mean_boot_t,
        Mean_boot_orig  = mean_boot_orig,
        Mean_boot_bench = mean_boot_bench,
        error_sd_boot   = err_sd,
        ran_eff_sd_boot = ran_sd
      )
    }

    # Skip the full-refit machinery below.
    return(.mse_megb_collate(
      boots_models = boots_models, tau_star = tau_star, pop_data = pop_data,
      dom_name = dom_name, corrected_bt = corrected_bt,
      benchmark_fn = benchmark_fn, domains = domains, error_sd = error_sd
    ))
  }

  if (bootstrap_refit == "leaves_only") {
    if (gbm_engine != "xgboost") {
      warning("bootstrap_refit='leaves_only' currently supports xgboost only; ",
              "falling back to 'full' for gbm_engine='", gbm_engine, "'.")
      # Falls through to the legacy full-refit block below. We DO NOT call
      # the .mse_megb_collate here; control returns to the wider function flow.
    } else {

      # Original booster.
      orig_booster <- model$boosting

      # Refresh updater nrounds: default to the model's original training
      # nrounds (refresh each tree in the ensemble once, in sequence).
      #
      # NOTE: xgboost's refresh updater is cascading - each tree's leaf
      # refresh uses residuals from the previously-refreshed earlier trees,
      # so leaf perturbations compound multiplicatively across the ensemble.
      # A diagnostic sweep (see diagnose_refresh_nrounds_sweep.R) showed
      # that the resulting MSE varies non-monotonically with nrounds: small
      # values under-estimate by failing to capture GB-refit variance
      # (the GB barely moves, the LMM trivially recovers the known u_d*
      # perturbation it was just handed), while large values over-estimate
      # via cascade artifact. The model's original nrounds is the most
      # defensible structural default - it matches what a naive reader of
      # the published REBB procedure (Messer & Schmid 2024) would expect -
      # but the variance estimate it produces may be inflated. Users
      # investigating uncertainty for sparse outcomes should consider this
      # an open methodological question and may wish to override via
      #   options(povmap.leaves_only.refresh_nrounds = N).
      .refresh_n_override <- getOption("povmap.leaves_only.refresh_nrounds")
      orig_nrounds <- if (!is.null(.refresh_n_override)) {
        as.integer(.refresh_n_override)
      } else {
        tryCatch(
          xgboost::xgb.get.num.boosted.rounds(orig_booster),
          error = function(e) {
            n <- gradient_params$nround
            if (is.null(n)) n <- gradient_params$nrounds
            if (is.null(n)) stop("Could not determine original nrounds for refresh.")
            n
          }
        )
      }

      # Feature matrices (numeric). cov_names was used by the original fit, so
      # column ordering will match.
      X_smp_mat <- as.matrix(smp_data[, cov_names, drop = FALSE])
      X_pop_mat <- as.matrix(pop_data[, cov_names, drop = FALSE])

      # LMM scaffolding (same as lmm_only branch).
      lmm_formula <- stats::as.formula(paste0("r ~ 1 + ", formula_random_effects))
      newdat_d    <- stats::setNames(
        data.frame(domains, stringsAsFactors = FALSE), dom_name
      )
      if (is.factor(smp_data[[dom_name]])) {
        newdat_d[[dom_name]] <- factor(domains, levels = levels(smp_data[[dom_name]]))
      }
      fit_data <- smp_data

      # collect_diag is defined at function scope at the top of mse_megb.
      if (collect_diag) {
        .D <- length(domains)
        gb_pop_d_mean_boot_mat <- matrix(NA_real_, nrow = .D, ncol = B)
        u_d_boot_mat           <- matrix(NA_real_, nrow = .D, ncol = B)
      }

      boots_models <- vector("list", B)

      for (i in seq_len(B)) {
        if (i %% max(1L, B %/% 10L) == 0L)
          message("leaves_only bootstrap iteration ", i, " of ", B)

        # 1) Refresh leaves on the bootstrap label. Tree structure unchanged.
        d_train <- xgboost::xgb.DMatrix(X_smp_mat, label = y_star_smp[, i])
        gb_refreshed <- tryCatch(
          suppressMessages(suppressWarnings(
            xgboost::xgb.train(
              params    = list(
                updater       = "refresh",
                process_type  = "update",
                refresh_leaf  = 1,
                objective     = "reg:squarederror"
              ),
              data      = d_train,
              nrounds   = orig_nrounds,
              xgb_model = orig_booster,
              verbose   = 0
            )
          )),
          error = function(e) {
            message("  refresh failed at iter ", i, ": ", conditionMessage(e))
            NULL
          }
        )
        if (is.null(gb_refreshed)) {
          boots_models[[i]] <- list(
            Mean_boot = NULL, Mean_boot_orig = NULL, Mean_boot_bench = NULL,
            error_sd_boot = NA_real_, ran_eff_sd_boot = NA_real_
          )
          next
        }

        # 2) Predict refreshed model on pop and smp.
        gb_pop_boot <- as.numeric(predict(gb_refreshed, X_pop_mat))
        gb_smp_boot <- as.numeric(predict(gb_refreshed, X_smp_mat))

        # 3) Refit LMM on residuals from refreshed predictions.
        fit_data$r <- y_star_smp[, i] - gb_smp_boot
        if (!is.null(smp_weights_vec)) fit_data$.megb_w <- smp_weights_vec
        lmer_boot <- tryCatch(
          suppressMessages(suppressWarnings(
            lme4::lmer(
              lmm_formula, data = fit_data, REML = TRUE,
              weights = if (!is.null(smp_weights_vec)) fit_data$.megb_w else NULL
            )
          )),
          error = function(e) NULL
        )
        if (is.null(lmer_boot)) {
          boots_models[[i]] <- list(
            Mean_boot = NULL, Mean_boot_orig = NULL, Mean_boot_bench = NULL,
            error_sd_boot = NA_real_, ran_eff_sd_boot = NA_real_
          )
          next
        }

        re_per_d <- stats::predict(lmer_boot, newdata = newdat_d,
                                   allow.new.levels = TRUE)
        fe_boot  <- as.numeric(lme4::fixef(lmer_boot))[1]
        u_d_boot <- as.numeric(re_per_d) - fe_boot




        # 4) Per-domain mean from refreshed pop predictions + BLUP.
        gb_pop_d_mean_boot <- as.numeric(tapply(gb_pop_boot, pop_dom_idx, mean))
        if (collect_diag) {
          gb_pop_d_mean_boot_mat[, i] <- gb_pop_d_mean_boot
          u_d_boot_mat[, i]           <- u_d_boot
        }
        mean_boot_t        <- gb_pop_d_mean_boot + u_d_boot
        mean_boot_orig     <- if (!is.null(corrected_bt)) corrected_bt(mean_boot_t)
        else mean_boot_t

        mean_boot_bench <- if (!is.null(benchmark_fn)) {
          bs_for_bench <- smp_data
          bs_for_bench$y_star <- y_star_smp[, i]
          benchmark_fn(mean_boot_orig, bs_for_bench)
        } else NULL

        ran_sd <- as.data.frame(lme4::VarCorr(lmer_boot))$sdcor[1]
        err_sd <- stats::sigma(lmer_boot)

        boots_models[[i]] <- list(
          Mean_boot       = mean_boot_t,
          Mean_boot_orig  = mean_boot_orig,
          Mean_boot_bench = mean_boot_bench,
          error_sd_boot   = err_sd,
          ran_eff_sd_boot = ran_sd
        )
      }

      if (collect_diag) {
        diag_env <- getOption("povmap.leaves_only.diag_env")
        assign("gb_pop_d_mean_boot_mat", gb_pop_d_mean_boot_mat, envir = diag_env)
        assign("u_d_boot_mat",           u_d_boot_mat,           envir = diag_env)
        assign("domains",                domains,                envir = diag_env)
      }

      return(.mse_megb_collate(
        boots_models = boots_models, tau_star = tau_star, pop_data = pop_data,
        dom_name = dom_name, corrected_bt = corrected_bt,
        benchmark_fn = benchmark_fn, domains = domains, error_sd = error_sd
      ))
    }
  }

  # ── Else: full refit (legacy path, both GB and LMM each iteration) ──────────
  my_estim_f <- function(x) {
    model_boot <- em_gb_lmm(
      Y                      = x$y_star,
      X                      = x[, .cov_names, drop = FALSE],
      dom_name               = .dom_name,
      data                   = x,
      pop_data               = .pop_data,
      seed                   = 1,
      gradient_params        = .gradient_params,
      formula_random_effects = .formula_re,
      initial_random_effects = 0,
      max_iterations         = 10,
      error_tolerance        = 1e-04,
      cov_names              = .cov_names,
      gbm_engine             = .gbm_engine,
      weights                = .smp_weights_vec
    )

    unit_level_predictions <- gbm_predict(
      model      = model_boot,
      smp_data   = .smp_data,
      pop_data   = .pop_data,
      Y          = .Y,
      dom_name   = .dom_name,
      gbm_engine = .gbm_engine,
      cov_names  = .cov_names
    )

    mean_preds <- unit_level_predictions$unit_pred_pop |>
      dplyr::group_by(dom_name) |>
      dplyr::summarise(Mean = mean(unit_preds)) |>
      as.data.frame()

    mean_boot_t <- mean_preds[, "Mean"]

    mean_boot_orig  <- if (!is.null(.corrected_bt)) .corrected_bt(mean_boot_t)
                       else mean_boot_t
    mean_boot_bench <- if (!is.null(.benchmark_fn)) .benchmark_fn(mean_boot_orig, x)
                       else NULL

    list(
      Mean_boot       = mean_boot_t,
      Mean_boot_orig  = mean_boot_orig,
      Mean_boot_bench = mean_boot_bench,
      error_sd_boot   = model_boot$error_sd,
      ran_eff_sd_boot = model_boot$ran_eff_sd
    )
  }
  # Lean environment: prevents furrr from serialising the full mse_megb frame to
  # each worker. Critically, excludes model$boosting (externalptr, non-exportable)
  # and large intermediate matrices (e_ij, u_i, pred_mat, boots_sample, etc.).
  environment(my_estim_f) <- list2env(
    list(
      .cov_names       = cov_names,
      .dom_name        = dom_name,
      .pop_data        = pop_data,
      .gradient_params = gradient_params,
      .formula_re      = formula_random_effects,
      .gbm_engine      = gbm_engine,
      .smp_data        = smp_data,
      .Y               = Y,
      .corrected_bt    = corrected_bt,
      .benchmark_fn    = benchmark_fn,
      .smp_weights_vec = smp_weights_vec
    ),
    parent = getNamespace("povmap")
  )

  if (is.null(bootstrap_cores) || bootstrap_cores <= 1) {
    boots_models <- vector("list", B)
    for (i in seq_len(B)) {
      message("Bootstrap iteration ", i, " of ", B)
      boots_models[[i]] <- my_estim_f(boots_sample[[i]])
    }
  } else {
    os_type <- Sys.info()[["sysname"]]
    if (os_type == "Windows") {
      future::plan(future::multisession, workers = bootstrap_cores)
    } else {
      future::plan(future::multicore, workers = bootstrap_cores)
    }
    old_max <- getOption("future.globals.maxSize")
    options(future.globals.maxSize = 4 * 1024^3)  # 4 GiB: pop_data can be large
    tryCatch(
      boots_models <- furrr::future_map(
        boots_sample, my_estim_f,
        .options  = furrr::furrr_options(seed = seed),
        .progress = TRUE
      ),
      finally = {
        options(future.globals.maxSize = old_max)
        future::plan(future::sequential)
      }
    )
  }

  .mse_megb_collate(
    boots_models = boots_models, tau_star = tau_star, pop_data = pop_data,
    dom_name = dom_name, corrected_bt = corrected_bt,
    benchmark_fn = benchmark_fn, domains = domains, error_sd = error_sd
  )
}

# Internal: collate per-iteration bootstrap outputs into the matrices and
# summary frames mse_megb returns. Used by both the lmm_only and full paths.
.mse_megb_collate <- function(boots_models, tau_star, pop_data, dom_name,
                              corrected_bt, benchmark_fn, domains, error_sd) {

  # Robust column-bind of per-iteration vectors of length D. Drops iterations
  # whose slot is NULL or whose length disagrees with D (with a warning).
  D <- nrow(tau_star)
  collate_boot <- function(slot) {
    cols <- lapply(boots_models, function(b) {
      v <- b[[slot]]
      if (is.null(v)) return(NULL)
      v <- as.numeric(v)
      if (length(v) != D) return(NULL)
      v
    })
    keep <- which(!vapply(cols, is.null, logical(1)))
    n_drop <- length(cols) - length(keep)
    if (n_drop > 0L)
      warning("Dropped ", n_drop, " of ", length(cols),
              " bootstrap iterations with NULL or wrong-length '", slot, "'.")
    if (length(keep) == 0L) return(list(mat = NULL, keep = integer(0)))
    list(mat = do.call(cbind, cols[keep]), keep = keep)
  }

  cb_t    <- collate_boot("Mean_boot")
  cb_orig <- collate_boot("Mean_boot_orig")
  tau_b_t      <- cb_t$mat
  tau_b_orig   <- cb_orig$mat
  keep_orig    <- cb_orig$keep
  if (is.null(tau_b_t))
    stop("All bootstrap iterations failed to return a valid Mean_boot vector.")

  MSE_estimates <- rowMeans((tau_star - tau_b_t)^2, na.rm = TRUE)
  MSE_estimates <- data.frame(unique(pop_data[dom_name]), Mean = MSE_estimates)
  rownames(MSE_estimates) <- NULL

  tau_star_orig       <- NULL
  tau_star_orig_unb   <- NULL
  tau_b_bench         <- NULL
  tau_star_orig_bench <- NULL
  MSE_bench_estimates <- NULL
  if (!is.null(corrected_bt)) {
    tau_star_orig     <- apply(tau_star, 2, corrected_bt)
    tau_star_orig_unb <- tau_star_orig[, keep_orig, drop = FALSE]
  }
  if (!is.null(corrected_bt) && !is.null(benchmark_fn)) {
    cb_bench    <- collate_boot("Mean_boot_bench")
    tau_b_bench <- cb_bench$mat
    if (!is.null(tau_b_bench) && !is.null(tau_star_orig)) {
      tau_star_orig_bench <- tau_star_orig[, cb_bench$keep, drop = FALSE]
      mse_bench_vec <- rowMeans((tau_star_orig_bench - tau_b_bench)^2, na.rm = TRUE)
      MSE_bench_estimates <- data.frame(unique(pop_data[dom_name]), Mean = mse_bench_vec)
      rownames(MSE_bench_estimates) <- NULL
    }
  }

  boot_error_sd        <- sapply(boots_models, getElement, "error_sd_boot")   |> unlist()
  boot_ran_eff_sd_boot <- sapply(boots_models, getElement, "ran_eff_sd_boot") |> unlist()

  list(
    call                 = match.call(),
    MSE_estimates        = MSE_estimates,
    MSE_bench_estimates  = MSE_bench_estimates,
    tau_b_t              = tau_b_t,
    tau_b_orig           = tau_b_orig,
    tau_b_bench          = tau_b_bench,
    tau_star             = tau_star,
    tau_star_orig        = tau_star_orig,
    tau_star_orig_unb    = tau_star_orig_unb,
    tau_star_orig_bench  = tau_star_orig_bench,
    domains              = domains,
    boot_error_sd        = boot_error_sd,
    boot_ran_eff_sd_boot = boot_ran_eff_sd_boot,
    error_sd_input       = error_sd
  )
}
