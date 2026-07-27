#' Mixed Effects Gradient Boosting for domain-level averages
#'
#' The function \code{megb} combines gradient boosting with linear mixed models
#' (via an EM algorithm) to estimate domain-level averages for small area
#' estimation. It wraps the \code{megb} function from the \pkg{MEGB} package,
#' returning a \code{povmap} object compatible with \code{write.excel},
#' \code{estimators}, and \code{map_plot}.
#'
#' Unlike \code{\link{xgb}}, MEGB operates directly at the unit level and
#' models area-level random effects explicitly through a linear mixed model.
#' No sub-domain identifier is required: the EM algorithm handles the
#' sample-to-population matching internally.
#'
#' @param fixed a two-sided linear formula describing the model. The dependent
#'   variable is on the left of \code{~} and the predictors on the right,
#'   separated by \code{+}. The domain identifier named in \code{domains} may
#'   appear in the formula (it is automatically excluded from the GB predictor
#'   matrix and used only as the LMM grouping variable).
#' @param smp_data a data frame containing all variables in \code{fixed} plus
#'   the domain identifier.
#' @param smp_weights character string naming the survey-weight variable in
#'   \code{smp_data}. Used only for domain-level aggregation of sample
#'   predictions (not passed to the GB engine). Defaults to \code{NULL}.
#' @param pop_data a data frame containing the domain identifier and all
#'   predictor variables from \code{fixed}.
#' @param pop_weights character string naming a population-count or expansion-
#'   weight variable in \code{pop_data}. When supplied, domain means are
#'   computed as weighted averages of unit-level predictions. Defaults to
#'   \code{NULL} (simple mean).
#' @param domains character string naming the domain variable in both
#'   \code{smp_data} and \code{pop_data}.
#' @param transformation transformation applied to the dependent variable before
#'   fitting. One of \code{"no"} (default), \code{"log"}, \code{"log.shift"},
#'   \code{"arcsin"}, \code{"logistic"}, or \code{"poisson"}. Unit-level
#'   population predictions are back-transformed before domain aggregation.
#' @param mse logical. If \code{TRUE} (default), runs the parametric bootstrap
#'   to quantify uncertainty. The form of the resulting summary is controlled
#'   by \code{mse_type}.
#'   The EM loop runs to convergence (relative log-likelihood change <
#'   1e-04) or for at most 25 iterations, whichever comes first. If
#'   \eqn{\sigma_u} collapses toward zero, that typically reflects
#'   identifiability between area-level covariates (which are constant
#'   within an area, and so collinear with the area indicator) and the
#'   area random effect, rather than a failure of EM. Diagnose by
#'   comparing specifications with and without the area-level covariates
#'   and by leave-areas-out cross-validation, rather than by truncating
#'   the EM loop.
#' @param bench_target one of \code{"fixed"} (default) or \code{"random"}.
#'   Controls how the benchmark target is treated *inside the bootstrap*. "fixed"
#'   rakes every replicate to the observed (province/level) target, so the
#'   benchmarked MSE/CI is the spread of the estimate around the realised truth at
#'   the observed target — the correct choice for benchmarking to an observed survey
#'   target. "random" recomputes the target from each replicate's own resampled
#'   survey, which makes the benchmarked estimate track its own bootstrap mean and
#'   collapses the level variance, giving too-narrow intervals (census-validated:
#'   well-sampled-area coverage 0.74 under "random" vs 0.96 under "fixed"). Only
#'   matters when \code{benchmark} is character (internal, survey-derived
#'   benchmarking); for numeric or data.frame benchmarks the target is a
#'   constant by construction.
#'   \code{"random"} (default): recompute the target on every bootstrap
#'   iteration from the bootstrap-simulated survey data. Both the target and
#'   the model estimate vary together, preserving the joint sampling
#'   correlation between them — which is the right thing for re-sampling-based
#'   inference about the unconditional variance of the benchmarked estimator.
#'   In practice this is also typically *narrower* than \code{"fixed"} for
#'   internal benchmarks: when survey resamples shift the target up they also
#'   shift the model estimate up, so the benchmark-adjustment shift \eqn{c}
#'   has less work to do, and \eqn{var(c) \propto var(T - WM)} stays small.
#'   \code{"fixed"}: hold the target at the value computed from the *observed*
#'   survey and reuse it on every iteration. The resulting variance is
#'   conditional on the target — \eqn{Var(\hat\theta_{bench} \mid target)}.
#'   Use this when the benchmark target really is a known constraint that you
#'   don't want to bootstrap over. For survey-derived (internal) benchmarks
#'   the target *is* itself a sample statistic, so this conditioning is a
#'   somewhat artificial inferential target — \code{"random"} is usually
#'   preferable.
#' @param bootstrap_refit one of \code{"leaves_only"} (default),
#'   \code{"full"}, or \code{"lmm_only"}.
#'   \code{"leaves_only"} (xgboost only) freezes the tree structure of the
#'   original GB across all bootstrap iterations and only refreshes leaf-node
#'   values via xgboost's \code{refresh} updater. This captures within-model
#'   leaf-prediction sampling uncertainty plus the LMM random-effect
#'   uncertainty, without exposing the SE to GB tree-selection variance —
#'   which empirically overstates uncertainty for sparse or boundary-near
#'   indicators (e.g. ownership rates, livestock counts, access indicators
#'   with many near-zero territories) where bootstrap-perturbed residuals
#'   drive fresh GB fits to substantially different trees. For dense
#'   indicators (continuous welfare, well-spread proportions) \code{"leaves_only"}
#'   gives results within ~20\% of \code{"full"}; for sparse indicators the
#'   two can differ by an order of magnitude. Recommended default and aligned
#'   with how most SAE methods report SEs (conditional on chosen model
#'   structure).
#'   \code{"full"}: refit both the gradient booster and the linear mixed model
#'   in every bootstrap iteration. This additionally captures GB
#'   tree-selection variance and is appropriate when that is the inferential
#'   target (e.g. inference about a model \emph{class} rather than a fitted
#'   model). For sparse indicators this can substantially overstate true
#'   predictive uncertainty because the residual bootstrap drives fresh GBs
#'   into pathological refits.
#'   \code{"lmm_only"}: refit only the linear mixed model on bootstrap residuals
#'   each iteration, treating the gradient-booster fit as fixed. This is the
#'   classical EBLUP parametric bootstrap (Prasad–Rao / Hall–Maiti) and yields
#'   the textbook leading-order MSE \eqn{g_1 = \gamma_d \sigma_e^2 / n_d}.
#'   It runs faster than the other modes but assumes the FE estimator is
#'   parametric and contributes negligible variance — an assumption that
#'   holds for low-dimensional linear FE models but \emph{not} for boosted
#'   trees, where it typically understates SE by an order of magnitude. Use
#'   only if you explicitly want the classical-EBLUP variance decomposition.
#' @param mse_type one of \code{"var"} (default) or \code{"mse"}. \code{"var"}
#'   reports the empirical bootstrap variance of back-transformed domain means
#'   and uses bootstrap quantiles (recentred at the point estimate) for the CI,
#'   matching \code{\link{xgb}}. Benchmarked CIs use the centred residual
#'   between benchmarked-bootstrap and bootstrap-truth, again as in \code{xgb}.
#'   \code{"mse"} reports the prediction MSE \eqn{E[(\hat\theta-\theta)^2]} —
#'   computed in transformed space and back-mapped via the delta method for the
#'   unbenchmarked estimator, and directly in original space for the
#'   benchmarked estimator — and uses a normal-approximation CI. \code{"mse"}
#'   is generally wider than \code{"var"} because it includes the irreducible
#'   variability of the random effect.
#' @param B number of parametric bootstrap iterations for MSE estimation.
#'   Defaults to \code{100}.
#' @param bootstrap_cores number of cores for parallel bootstrap. \code{0} or
#'   \code{1} runs sequentially. Defaults to \code{0}.
#' @param conf_level confidence level for the confidence interval. Defaults to
#'   \code{0.95}.
#' @param gradient_params named list of gradient boosting hyperparameters passed
#'   to the chosen \code{gbm_engine}. When \code{NULL}, engine defaults are
#'   used. For XGBoost the relevant names are \code{eta}, \code{nrounds},
#'   \code{max_depth}, \code{subsample}, etc.
#' @param gbm_engine gradient boosting engine: \code{"xgboost"} (default),
#'   \code{"lightgbm"}, or \code{"catboost"}.
#' @param na.rm if \code{TRUE}, rows with \code{NA} values are removed from
#'   both \code{smp_data} and \code{pop_data} before fitting. Defaults to
#'   \code{FALSE}.
#' @param seed integer seed for reproducibility. Defaults to \code{123}.
#' @param ... additional arguments forwarded to \code{MEGB::megb}.
#'
#' @return An object of class \code{c("megb", "xgb", "povmap")} with elements:
#'   \describe{
#'     \item{\code{ind}}{data frame of domain-level point estimates
#'       (\code{Domain}, \code{Mean}).}
#'     \item{\code{var}}{data frame of delta-method-corrected MSE estimates
#'       (\code{Domain}, \code{Mean}), or \code{NULL} if \code{mse = FALSE}.}
#'     \item{\code{CI}}{data frame of normal-approximation confidence intervals
#'       (\code{Domain}, \code{Lower}, \code{Upper}), or \code{NULL} if
#'       \code{mse = FALSE}.}
#'     \item{\code{yhat}}{data frame of back-transformed sample-unit predictions
#'       (\code{obs_id}, \code{hat}).}
#'     \item{\code{model}}{the fitted GB model object (xgboost booster when
#'       \code{gbm_engine = "xgboost"}).}
#'     \item{\code{megb_model}}{the full MEGB model list including both the
#'       booster and the LME4 random-effects model.}
#'     \item{\code{gbm_engine}}{the engine used.}
#'     \item{\code{smp_data}}{the original sample data.}
#'     \item{\code{out_call}}{the matched call.}
#'     \item{\code{transformation}}{the transformation string.}
#'     \item{\code{framework}}{internal framework list with domain metadata.}
#'   }
#'
#' @seealso \code{\link{xgb}}, \code{\link{write.excel}},
#'   \code{\link{estimators}}
#'
#' @references
#' Merfeld, J. D., Dang, H., & Newhouse, D. (2025). Improving Estimates of
#' Mean Welfare and Uncertainty in Developing Countries. The World Bank.
#'
#' @importFrom stats qnorm tapply
#' @export
#'
#' @examples
#' \donttest{
#' data("eusilcA_pop")
#' data("eusilcA_smp")
#'
#' megb_model <- megb(
#'   fixed   = eqIncome ~ eqsize + cash + self_empl + unempl_ben +
#'             age_ben + surv_ben + sick_ben + dis_ben + rent +
#'             fam_allow + house_allow + cap_inv + tax_adj,
#'   smp_data = eusilcA_smp,
#'   pop_data = eusilcA_pop,
#'   domains  = "district",
#'   gradient_params = list(eta = 0.1, nrounds = 100, max_depth = 3)
#' )
#'
#' estimators(megb_model, indicator = "Mean", MSE = TRUE, CV = TRUE)
#' }

megb <- function(fixed,
                 smp_data,
                 smp_weights      = NULL,
                 pop_data,
                 pop_weights      = NULL,
                 domains,
                 transformation   = "no",
                 mse              = TRUE,
                 mse_type         = c("var", "mse"),
                 bootstrap_refit  = c("leaves_only", "full", "lmm_only"),
                 bench_target     = c("fixed", "random"),
                 B                = 100,
                 bootstrap_cores  = 0,
                 conf_level       = 0.95,
                 gradient_params  = NULL,
                 gbm_engine       = "xgboost",
                 na.rm            = FALSE,
                 seed             = 123,
                 benchmark        = NULL,
                 benchmark_type   = "ratio",
                 benchmark_level  = NULL,
                 benchmark_weights = NULL,
                 weightedBS       = TRUE,
                 cpus             = NULL,
                 ...) {
  # weightedBS: when TRUE and smp_weights is supplied, the bootstrap residual
  # sampling uses probabilities proportional to weights (mirroring xgb's
  # weightedBS semantics). When FALSE, residuals are sampled with equal
  # probability (the historical default before weights were plumbed through).

  out_call        <- match.call()
  mse_type        <- match.arg(mse_type)
  bootstrap_refit <- match.arg(bootstrap_refit)
  bench_target    <- match.arg(bench_target)

  # Local %||% so the diagnostic message below isn't fragile to NULLs.
  `%||%` <- function(a, b) if (is.null(a)) b else a

  if (!is.null(cpus)) bootstrap_cores <- cpus

  if (is.null(benchmark_weights) && !is.null(smp_weights))
    benchmark_weights <- smp_weights

  # ── 1. Framework ─────────────────────────────────────────────────────────────
  fwk <- framework_megb(
    fixed             = fixed,
    smp_data          = smp_data,
    pop_data          = pop_data,
    smp_weights       = smp_weights,
    pop_weights       = pop_weights,
    domains           = domains,
    transformation    = transformation,
    conf_level        = conf_level,
    na.rm             = na.rm,
    benchmark_weights = benchmark_weights
  )

  # ── 1b. Local benchmark helper ────────────────────────────────────────────────
  # Route through benchmark_xgb_level for level-based benchmarking so megb
  # supports the same benchmark_type set as xgb (notably "logit_raking",
  # which benchmark_ebp_level does not implement and would silently return
  # all-NA Mean_bench).
  add_benchmark_megb <- function(x, benchmark_level, fwk, fixed,
                                 benchmark, benchmark_type, domain_labels = NULL) {
    # Carry domain labels so benchmark_xgb_level can assert positional alignment
    # with unique(pop_data[[domains]]) (guards the key-alignment fix above).
    point_estim       <- if (!is.null(domain_labels))
                           list(ind = data.frame(Mean = x, Domain = domain_labels,
                                                 stringsAsFactors = FALSE))
                         else list(ind = data.frame(Mean = x))
    if (is.null(benchmark_level)) {
      point_estim$ind <- benchmark_ebp_national(
        point_estim    = point_estim,
        framework      = fwk,
        fixed          = fixed,
        benchmark      = benchmark,
        benchmark_type = benchmark_type)
    } else {
      point_estim$ind <- benchmark_xgb_level(
        point_estim     = point_estim,
        framework       = fwk,
        fixed           = fixed,
        benchmark       = benchmark,
        benchmark_type  = benchmark_type,
        benchmark_level = benchmark_level)
    }
    point_estim$ind$Mean_bench
  }
  # Lean env: only package functions needed, no local-scope capture
  environment(add_benchmark_megb) <- getNamespace("povmap")

  # ── 2. Transformation functions ───────────────────────────────────────────────
  if (transformation == "arcsin") {
    transform_outcome      <- arcsin_transform
    back_transform_outcome <- arcsin_transform_back
  } else if (transformation == "log") {
    transform_outcome      <- log_transform
    back_transform_outcome <- log_transform_back
  } else if (transformation == "log.shift") {
    ls_lambda              <- if (min(fwk$Y_smp) <= 0) abs(min(fwk$Y_smp)) + 1 else 0
    transform_outcome      <- function(y) list(y = log(y + ls_lambda), shift = NULL)
    back_transform_outcome <- function(y) exp(y) - ls_lambda
  } else if (transformation == "logistic") {
    transform_outcome      <- logit_transform_epsilon
    back_transform_outcome <- logit_transform_back
  } else if (transformation == "poisson") {
    if (min(fwk$Y_smp) < 0)
      stop("Outcome must be non-negative when using the poisson transformation")
    transform_outcome      <- poisson_transform
    back_transform_outcome <- poisson_transform_back
  } else {
    transform_outcome      <- no_transform
    back_transform_outcome <- no_transform_back
  }

  # ── 3. Transform outcome ──────────────────────────────────────────────────────
  Y_smp_t <- transform_outcome(fwk$Y_smp)$y

  # ── 4. Prepare data for MEGB::megb ───────────────────────────────────────────
  # Domain must be a factor; align factor levels across sample and population
  all_levels <- union(
    as.character(unique(fwk$smp_data[[domains]])),
    as.character(unique(fwk$pop_data[[domains]]))
  )

  smp_megb           <- fwk$smp_data[, c(domains, fwk$covariates), drop = FALSE]
  smp_megb[[domains]] <- factor(smp_megb[[domains]], levels = all_levels)

  pop_megb           <- fwk$pop_data[, c(domains, fwk$covariates), drop = FALSE]
  pop_megb[[domains]] <- factor(pop_megb[[domains]], levels = all_levels)

  X_smp <- smp_megb[, fwk$covariates, drop = FALSE]

  # ── 5. Fit MEGB model (MSE computed separately below after corrected_bt is known)
  megb_fit <- megb_em(
    Y               = Y_smp_t,
    X               = X_smp,
    dom_name        = domains,
    smp_data        = smp_megb,
    pop_data        = pop_megb,
    gradient_params = gradient_params,
    na.rm           = FALSE,
    seed            = seed,
    mse             = FALSE,
    gbm_engine      = gbm_engine,
    smp_weights_vec = fwk$smp_weights_vec,
    ...
  )

  # ── 5b. Jensen-corrected back-transformation ─────────────────────────────────
  # E[g^{-1}(hat_t + e)] > g^{-1}(hat_t) for convex g^{-1}. We correct using
  # the MEGB-estimated error variance sigma2_e. Analytical formulae follow
  # expected_untransformed_mean() in point_estimation.R (9th-order Taylor for
  # arcsin; exact lognormal formula for log/log.shift/poisson).
  sigma2_e <- tryCatch({
    err_sd <- megb_fit$megb_model$error_sd
    if (!is.null(err_sd) && !is.na(err_sd) && err_sd > 0) err_sd^2 else 0
  }, error = function(e) 0)

  # corrected_bt uses a lean captured environment (only scalar values + package ns)
  # so parallel workers don't get the full megb() scope serialised (~1.5 GiB).
  corrected_bt <- function(hat_t) {
    if (.sigma2 == 0 || .transf %in% c("no", "logistic")) {
      return(switch(.transf,
        "no"        = hat_t,
        "log"       = exp(hat_t),
        "log.shift" = exp(hat_t) - .ls_lam,
        "arcsin"    = arcsin_transform_back(hat_t),
        "logistic"  = logit_transform_back(hat_t),
        "poisson"   = expm1(hat_t),
        hat_t
      ))
    }
    switch(.transf,
      "log"       = exp(hat_t + 0.5 * .sigma2),
      "log.shift" = exp(hat_t + 0.5 * .sigma2) - .ls_lam,
      "poisson"   = pmax(0, expm1(hat_t + 0.5 * .sigma2)),
      "arcsin"    = {
        term1 <- arcsin_transform_back(hat_t)
        dy2dx <- -2 * sin(hat_t)^2 + 2 * cos(hat_t)^2
        term1 +
          0.5       * dy2dx * .sigma2 -
          (1/24)    * (-4  * dy2dx) * 3   * .sigma2^2 +
          (1/720)   * (-16 * dy2dx)       * .sigma2^3 +
          (1/40320) * (-64 * dy2dx)       * .sigma2^4
      },
      hat_t
    )
  }
  environment(corrected_bt) <- list2env(
    list(
      .sigma2 = sigma2_e,
      .transf = transformation,
      .ls_lam = if (transformation == "log.shift") ls_lambda else NULL
    ),
    parent = getNamespace("povmap")
  )

  # ── 6. Back-transform unit predictions and aggregate to domain level ──────────
  unit_preds_pop      <- megb_fit$unit_preds_all       # columns: unit_preds, dom_name
  unit_preds_pop$hat  <- corrected_bt(unit_preds_pop$unit_preds)

  if (!is.null(pop_weights)) {
    pw           <- fwk$pop_weights_vec                # same row order as pop_megb
    domain_means <- tapply(unit_preds_pop$hat * pw, unit_preds_pop$dom_name, sum) /
                    tapply(pw,                   unit_preds_pop$dom_name, sum)
  } else {
    domain_means <- tapply(unit_preds_pop$hat, unit_preds_pop$dom_name, mean)
  }

  ind <- data.frame(
    Domain = names(domain_means),
    Mean   = as.numeric(domain_means),
    stringsAsFactors = FALSE
  )

  # ── Benchmark key-alignment fix ───────────────────────────────────────────────
  # benchmark_xgb_level() positionally zips point_estim$ind$Mean against
  # unique(pop_data[[domains]]) (and the per-domain population weights). megb builds
  # `ind` via tapply(), which orders domains ALPHABETICALLY; pop_data is in
  # first-appearance order. When the two differ, every misplaced territory is raked
  # under the wrong benchmark-level target, silently corrupting the benchmarked
  # estimates (observed: ~90% of territories mis-assigned, with non-monotone
  # within-level redistribution and spurious zeros). xgb is unaffected because it
  # builds `ind` in pop_data order natively. Reorder `ind` to pop_data's domain
  # order so the positional contract holds; the assertion in benchmark_xgb_level
  # turns any residual mismatch into a loud error rather than silent corruption.
  .pop_dom_order <- as.character(unique(fwk$pop_data[[domains]]))
  ind <- ind[match(.pop_dom_order, as.character(ind$Domain)), , drop = FALSE]
  rownames(ind) <- NULL

  # ── 6b. Benchmarking of point estimates ───────────────────────────────────────
  ind_bench <- NULL
  if (!is.null(benchmark)) {
    bench_vals <- add_benchmark_megb(
      x               = ind$Mean,
      domain_labels   = ind$Domain,
      benchmark_level = benchmark_level,
      fwk             = fwk,
      fixed           = fixed,
      benchmark       = benchmark,
      benchmark_type  = benchmark_type
    )
    # Guard against silent benchmark failures (unsupported benchmark_type,
    # missing benchmark_level column, etc.) which previously surfaced as
    # all-empty Mean_bench / Var_bench / Lower_bench columns in write.excel.
    if (is.null(bench_vals) || length(bench_vals) != length(ind$Mean) ||
        all(is.na(bench_vals))) {
      stop(sprintf(
        "Benchmarking returned %s. Common causes: benchmark_type '%s' is not supported by the level-benchmarker (supported: 'ratio', 'ratio_bound', 'ratio_complement', 'logit_raking'); benchmark_level '%s' missing from smp_data or pop_data; or fixed[2] '%s' missing from smp_data.",
        if (is.null(bench_vals)) "NULL"
          else if (all(is.na(bench_vals))) "all-NA"
          else paste0("length ", length(bench_vals),
                      " (expected ", length(ind$Mean), ")"),
        as.character(benchmark_type),
        as.character(benchmark_level),
        as.character(fixed[[2]])
      ))
    }
    ind_bench <- data.frame(
      Domain = ind$Domain,
      Mean   = bench_vals,
      stringsAsFactors = FALSE
    )
  }

  # ── 6c. Build benchmark_fn closure for the bootstrap ─────────────────────────
  # When bench_target = "fixed" AND benchmark is character (internal), we
  # precompute the target ONCE from the observed survey using the same
  # survey-weighted-mean formula benchmark_xgb_level applies for internal
  # benchmarking. The closure then uses this fixed value instead of
  # recomputing from boot_smp_data$y_star each iteration, yielding a
  # conditional bench variance Var(theta_hat_bench | target = observed)
  # rather than the unconditional one that double-counts survey sampling
  # uncertainty on the target side.
  fixed_bm_val <- NULL
  if (!is.null(benchmark) && mse && bench_target == "fixed" &&
      !is.numeric(benchmark) && !is.data.frame(benchmark)) {
    bench_w_var <- if (!is.null(fwk$benchmark_weights)) fwk$benchmark_weights
                   else fwk$smp_weights
    if (!is.null(benchmark_level)) {
      grp_obs   <- fwk$smp_data[[benchmark_level]]
      y_obs     <- fwk$Y_smp
      w_obs     <- if (!is.null(bench_w_var)) fwk$smp_data[[bench_w_var]]
                   else rep(1, length(y_obs))
      lev_means <- tapply(seq_along(y_obs), grp_obs, function(idx)
        sum(y_obs[idx] * w_obs[idx]) / sum(w_obs[idx]))
      df <- data.frame(as.character(names(lev_means)),
                       as.numeric(lev_means),
                       stringsAsFactors = FALSE)
      names(df) <- c(benchmark_level, "Mean")
      fixed_bm_val <- df
    } else {
      # National-level internal benchmarking: scalar target.
      w_obs <- if (!is.null(bench_w_var)) fwk$smp_data[[bench_w_var]]
               else rep(1, length(fwk$Y_smp))
      fixed_bm_val <- c(Mean = sum(fwk$Y_smp * w_obs) / sum(w_obs))
    }
  }

  if (!is.null(benchmark) && mse) {
    # Use a lean captured environment: explicit values only, parent = package ns.
    # This prevents serialising the full megb() scope to each parallel worker.
    benchmark_fn <- function(est_orig, boot_smp_data) {
      # Branch order: explicit user-supplied target → precomputed fixed target
      # (internal + bench_target="fixed") → recompute from bootstrap survey
      # (random target).
      # Random recompute MUST use the same survey-weighted mean as the point
      # benchmark (benchmark_xgb_level uses w = framework$smp_weights). An
      # unweighted mean would produce a different statistic and jitter far
      # more than the weighted target under skewed survey weights, widening
      # bench CIs spuriously.
      bm_val <- if (is.numeric(.bm) || is.data.frame(.bm)) {
        .bm
      } else if (!is.null(.fixed_bm)) {
        .fixed_bm
      } else if (is.null(.bm_level)) {
        w_obs <- if (!is.null(.bm_w_var)) boot_smp_data[[.bm_w_var]]
                 else rep(1, nrow(boot_smp_data))
        y_bt  <- .corr_bt(boot_smp_data$y_star)
        c(Mean = sum(y_bt * w_obs) / sum(w_obs))
      } else {
        y_bt  <- .corr_bt(boot_smp_data$y_star)
        grp   <- boot_smp_data[[.bm_level]]
        w_obs <- if (!is.null(.bm_w_var)) boot_smp_data[[.bm_w_var]]
                 else rep(1, length(y_bt))
        lev_means <- tapply(seq_along(y_bt), grp, function(idx)
          sum(y_bt[idx] * w_obs[idx]) / sum(w_obs[idx]))
        df <- data.frame(as.character(names(lev_means)),
                         as.numeric(lev_means),
                         stringsAsFactors = FALSE)
        names(df) <- c(.bm_level, "Mean")
        df
      }
      # benchmark_xgb_level discards non-numeric benchmark inputs and instead
      # recomputes from framework$smp_data (because of an `if (!is.numeric)`
      # branch upstream). For a data.frame target like our precomputed
      # .fixed_bm or the per-iteration random target above, that means our
      # value would be silently ignored. Convert to a named numeric vector,
      # which goes through the branch that actually uses the input.
      if (is.data.frame(bm_val) && !is.null(.bm_level)) {
        bm_val <- stats::setNames(as.numeric(bm_val[["Mean"]]),
                                  as.character(bm_val[[.bm_level]]))
      }
      .add_bm(
        x               = est_orig,
        benchmark_level = .bm_level,
        fwk             = .fwk,
        fixed           = .fixed,
        benchmark       = bm_val,
        benchmark_type  = .bm_type
      )
    }
    environment(benchmark_fn) <- list2env(
      list(
        .bm       = benchmark,
        .bm_level = benchmark_level,
        .bm_type  = benchmark_type,
        .bm_w_var = if (!is.null(fwk$benchmark_weights)) fwk$benchmark_weights
                    else fwk$smp_weights,
        .fwk      = fwk,
        .fixed    = fixed,
        .corr_bt  = corrected_bt,
        .add_bm   = add_benchmark_megb,
        .fixed_bm = fixed_bm_val
      ),
      parent = getNamespace("povmap")
    )
  } else {
    benchmark_fn <- NULL
  }

  # ── 6d. Parametric bootstrap MSE (now that corrected_bt/benchmark_fn are ready)
  # pop_data_proc / smp_data are one-hot encoded and lack the benchmark_level
  # column. Re-attach it from fwk$pop_data and fwk$smp_data so the bootstrap
  # samples (now derived from smp_data, not a pop subsample) carry it through
  # to benchmark_fn, which calls tapply(y_bt, boot_smp_data[[bm_level]], ...).
  boot_pop_data <- megb_fit$pop_data_proc
  boot_smp_data <- megb_fit$inp_smp_data$smp_data
  if (!is.null(benchmark_level) && mse) {
    if (!benchmark_level %in% colnames(boot_pop_data)) {
      bm_pop <- fwk$pop_data[[benchmark_level]]
      if (!is.null(bm_pop) && length(bm_pop) == nrow(boot_pop_data))
        boot_pop_data[[benchmark_level]] <- bm_pop
      else
        warning("benchmark_level '", benchmark_level,
                "' not found in pop_data or row count mismatch — ",
                "benchmarked bootstrap MSE will be skipped.")
    }
    if (!benchmark_level %in% colnames(boot_smp_data)) {
      bm_smp <- fwk$smp_data[[benchmark_level]]
      if (!is.null(bm_smp) && length(bm_smp) == nrow(boot_smp_data))
        boot_smp_data[[benchmark_level]] <- bm_smp
      else
        warning("benchmark_level '", benchmark_level,
                "' not found in smp_data or row count mismatch — ",
                "benchmarked bootstrap MSE will be skipped.")
    }
  }

  mse_estimated <- NULL
  if (mse) {
    message("Bootstrap with ", B, " iterations has started")
    # Attach weights as a column on boot_smp_data to guarantee row alignment.
    # fwk$smp_weights_vec was built from fwk$smp_data, but boot_smp_data
    # (= megb_fit$inp_smp_data$smp_data) may have a different row count or
    # ordering after framework + megb_em processing. Pulling weights from a
    # column ensures they always match the smp_data rows downstream consumers
    # iterate over.
    if (!is.null(fwk$smp_weights_vec)) {
      if (length(fwk$smp_weights_vec) == nrow(boot_smp_data)) {
        boot_smp_data$.megb_w <- as.numeric(fwk$smp_weights_vec)
      } else {
        warning("smp_weights_vec length (", length(fwk$smp_weights_vec),
                ") does not match boot_smp_data row count (",
                nrow(boot_smp_data), "). Falling back to unweighted bootstrap.")
        boot_smp_data$.megb_w <- NULL
      }
    }
    mse_estimated <- mse_megb(
      Y                      = megb_fit$inp_smp_data$target_var,
      X                      = megb_fit$X_proc,
      dom_name               = domains,
      smp_data               = boot_smp_data,
      model                  = megb_fit$megb_model,
      error_sd               = megb_fit$megb_model$error_sd,
      pop_data               = boot_pop_data,
      B                      = B,
      initial_random_effects = 0,
      ErrorTolerance         = 0.0001,
      MaxIterations          = 10,
      cov_names              = megb_fit$cov_names_proc,
      gradient_params        = gradient_params,
      formula_random_effects = megb_fit$formula_random_effects,
      bootstrap_cores        = bootstrap_cores,
      seed                   = seed,
      gbm_engine             = gbm_engine,
      unit_pred_smp          = megb_fit$unit_pred_smp,
      gb_smp                 = megb_fit$gb_smp,
      unit_preds             = megb_fit$unit_preds_all,
      bootstrap_refit        = bootstrap_refit,
      corrected_bt           = corrected_bt,
      benchmark_fn           = benchmark_fn,
      smp_weights_col        = ".megb_w",
      weightedBS             = weightedBS,
      # Same per-cell population weights the point estimate aggregates with
      # (megb.R: domain_means = sum(hat*pw)/sum(pw) over pop cells). Passed so the
      # bootstrap's per-cell-then-aggregate domain means use identical weighting;
      # NULL reproduces the historical unweighted bootstrap aggregation.
      pop_weights_vec        = if (!is.null(pop_weights)) fwk$pop_weights_vec else NULL
    )
  }

  # ── 7. Variance / MSE and CI ─────────────────────────────────────────────────
  # mse_type = "var" (default, mirrors xgb): empirical bootstrap variance and
  #   quantile-based CIs in original space. Benchmarked CIs use the centred
  #   residual (boot_bench - boot_truth) recentred at the point benchmark, as
  #   in xgb.R.
  # mse_type = "mse": prediction MSE in original space. Unbenchmarked MSE is
  #   computed in transformed space from the bootstrap and back-mapped via the
  #   delta method; benchmarked MSE is computed directly in original space.
  #   CIs use a normal approximation around the point estimate.
  var_df       <- NULL
  ci_df        <- NULL
  var_bench_df <- NULL
  ci_bench_df  <- NULL

  if (mse && !is.null(mse_estimated)) {

    alpha <- 1 - conf_level
    z_val <- qnorm(1 - alpha / 2)

    # Helper: clamp CIs to the valid outcome support for bounded transformations.
    clamp_ci <- function(df) {
      if (transformation %in% c("arcsin", "logistic")) {
        df$Lower <- pmax(df$Lower, 0); df$Upper <- pmin(df$Upper, 1)
      } else if (transformation %in% c("poisson", "log", "log.shift")) {
        df$Lower <- pmax(df$Lower, 0)
      }
      df
    }

    boot_domains <- as.character(mse_estimated$domains)

    if (mse_type == "var") {

      # ── 7a-var. Unbenchmarked: centred residuals of (boot - boot_truth) ──
      # For a mixed-effects estimator the target of inference is the realised
      # domain mean (which depends on u_d), so the natural uncertainty summary
      # is the prediction error against the bootstrap truth, not the variance
      # of the bootstrap estimates alone. This is symmetric with how the
      # benchmarked CI is computed below and with the "mse" path.
      tau_b_orig <- mse_estimated$tau_b_orig         # D × B' (original space)
      truth_u    <- mse_estimated$tau_star_orig_unb  # truth restricted to same B'

      resid_u <- tau_b_orig - truth_u
      resid_u_centred <- resid_u - rowMeans(resid_u, na.rm = TRUE)

      var_boot <- apply(resid_u_centred, 1, var,      na.rm = TRUE)
      lo_boot  <- apply(resid_u_centred, 1, quantile, probs = alpha / 2,     na.rm = TRUE)
      hi_boot  <- apply(resid_u_centred, 1, quantile, probs = 1 - alpha / 2, na.rm = TRUE)

      ind$Domain <- as.character(ind$Domain)
      idx        <- match(ind$Domain, boot_domains)

      var_df <- data.frame(Domain = ind$Domain, Mean = var_boot[idx],
                           stringsAsFactors = FALSE)

      ci_df <- data.frame(
        Domain = ind$Domain,
        Lower  = ind$Mean + lo_boot[idx],
        Upper  = ind$Mean + hi_boot[idx],
        stringsAsFactors = FALSE
      )
      ci_df <- clamp_ci(ci_df)

      # ── 7b-var. Benchmarked: var/quantiles of (boot_bench - boot_truth) ──
      # Both matrices are restricted to the same kept iterations by mse_megb,
      # so column counts match without further index gymnastics.
      if (!is.null(benchmark) && !is.null(mse_estimated$tau_b_bench) &&
          !is.null(mse_estimated$tau_star_orig_bench)) {
        tau_b_bench <- mse_estimated$tau_b_bench
        truth_bench <- mse_estimated$tau_star_orig_bench
        if (ncol(tau_b_bench) > 0L) {
          resid_mat     <- tau_b_bench - truth_bench
          row_means     <- rowMeans(resid_mat, na.rm = TRUE)
          resid_centred <- resid_mat - row_means

          var_bench_vec <- apply(resid_centred, 1, var,      na.rm = TRUE)
          lo_bench_vec  <- apply(resid_centred, 1, quantile, probs = alpha / 2,     na.rm = TRUE)
          hi_bench_vec  <- apply(resid_centred, 1, quantile, probs = 1 - alpha / 2, na.rm = TRUE)

          ind_bench$Domain <- as.character(ind_bench$Domain)
          idx_b <- match(ind_bench$Domain, boot_domains)

          var_bench_df <- data.frame(Domain = ind_bench$Domain,
                                     Mean   = var_bench_vec[idx_b],
                                     stringsAsFactors = FALSE)
          ci_bench_df <- data.frame(
            Domain = ind_bench$Domain,
            Lower  = ind_bench$Mean + lo_bench_vec[idx_b],
            Upper  = ind_bench$Mean + hi_bench_vec[idx_b],
            stringsAsFactors = FALSE
          )
          ci_bench_df <- clamp_ci(ci_bench_df)
        }
      }

    } else {  # mse_type == "mse"

      # ── 7a-mse. Unbenchmarked prediction MSE (delta-method back-map) ──
      mse_raw              <- mse_estimated$MSE_estimates
      colnames(mse_raw)[1] <- "Domain"
      mse_raw$Domain       <- as.character(mse_raw$Domain)
      colnames(mse_raw)[2] <- "MSE_t"

      mean_t        <- as.data.frame(megb_fit$Indicators)
      colnames(mean_t)[colnames(mean_t) == "dom_name"] <- "Domain"
      colnames(mean_t)[colnames(mean_t) == "Mean"]     <- "Mean_t"
      mean_t$Domain <- as.character(mean_t$Domain)

      ind$Domain <- as.character(ind$Domain)
      merged     <- merge(ind,    mse_raw, by = "Domain")
      merged     <- merge(merged, mean_t,  by = "Domain")

      delta_factor <- switch(transformation,
        "no"        = rep(1, nrow(merged)),
        "log"       = exp(merged$Mean_t)^2,
        "log.shift" = exp(merged$Mean_t)^2,
        "arcsin"    = sin(2 * merged$Mean_t)^2,
        "logistic"  = { p <- back_transform_outcome(merged$Mean_t); (p * (1 - p))^2 },
        "poisson"   = exp(merged$Mean_t)^2
      )
      merged$MSE_orig <- merged$MSE_t * delta_factor

      var_df <- data.frame(Domain = merged$Domain, Mean = merged$MSE_orig,
                           stringsAsFactors = FALSE)
      ci_df <- data.frame(
        Domain = merged$Domain,
        Lower  = merged$Mean - z_val * sqrt(pmax(merged$MSE_orig, 0)),
        Upper  = merged$Mean + z_val * sqrt(pmax(merged$MSE_orig, 0)),
        stringsAsFactors = FALSE
      )
      ci_df <- clamp_ci(ci_df)
      ind   <- data.frame(Domain = merged$Domain, Mean = merged$Mean,
                          stringsAsFactors = FALSE)

      # ── 7b-mse. Benchmarked prediction MSE (already in original space) ──
      if (!is.null(benchmark) && !is.null(mse_estimated$MSE_bench_estimates)) {
        mse_bench_raw              <- mse_estimated$MSE_bench_estimates
        colnames(mse_bench_raw)[1] <- "Domain"
        mse_bench_raw$Domain       <- as.character(mse_bench_raw$Domain)
        colnames(mse_bench_raw)[2] <- "MSE_bench"

        ind_bench$Domain <- as.character(ind_bench$Domain)
        merged_b         <- merge(ind_bench, mse_bench_raw, by = "Domain")

        var_bench_df <- data.frame(Domain = merged_b$Domain, Mean = merged_b$MSE_bench,
                                   stringsAsFactors = FALSE)
        ci_bench_df <- data.frame(
          Domain = merged_b$Domain,
          Lower  = merged_b$Mean - z_val * sqrt(pmax(merged_b$MSE_bench, 0)),
          Upper  = merged_b$Mean + z_val * sqrt(pmax(merged_b$MSE_bench, 0)),
          stringsAsFactors = FALSE
        )
        ci_bench_df <- clamp_ci(ci_bench_df)
        ind_bench   <- data.frame(Domain = merged_b$Domain, Mean = merged_b$Mean,
                                  stringsAsFactors = FALSE)
      }
    }
  }

  # ── 8. Sample-level predictions (back-transformed) ───────────────────────────
  unit_pred_smp_bt <- corrected_bt(megb_fit$unit_pred_smp)
  yhat <- data.frame(
    obs_id = seq_along(unit_pred_smp_bt),
    hat    = as.numeric(unit_pred_smp_bt)
  )

  # Store in framework for summary use
  fwk$unit_pred_smp <- as.numeric(unit_pred_smp_bt)

  # ── 9. Model reference ────────────────────────────────────────────────────────
  # For xgboost: expose the booster so add_model_megb can compute importance
  if (gbm_engine == "xgboost") {
    model_ref <- megb_fit$megb_model$boosting
  } else {
    model_ref <- megb_fit$megb_model
  }

  # ── 9b. Align result structure with xgb so shared write.excel/estimators logic works ──
  # xgb stores Mean_bench as a column in $ind, Var_bench in $var, Lower/Upper_bench in $CI.
  if (!is.null(ind_bench)) {
    ind <- merge(
      ind,
      data.frame(Domain = ind_bench$Domain, Mean_bench = ind_bench$Mean,
                 stringsAsFactors = FALSE),
      by = "Domain", all.x = TRUE
    )
  }
  if (!is.null(var_df) && !is.null(var_bench_df)) {
    var_df <- merge(
      var_df,
      data.frame(Domain = var_bench_df$Domain, Var_bench = var_bench_df$Mean,
                 stringsAsFactors = FALSE),
      by = "Domain", all.x = TRUE
    )
  }
  if (!is.null(ci_df) && !is.null(ci_bench_df)) {
    ci_bench_renamed <- ci_bench_df
    names(ci_bench_renamed)[names(ci_bench_renamed) == "Lower"] <- "Lower_bench"
    names(ci_bench_renamed)[names(ci_bench_renamed) == "Upper"] <- "Upper_bench"
    ci_df <- merge(ci_df, ci_bench_renamed, by = "Domain", all.x = TRUE)
  }

  # ── 9c. Bootstrap diagnostics ────────────────────────────────────────────────
  # Surface boot_ran_eff_sd / boot_error_sd alongside their non-bootstrap
  # counterparts so the caller can diagnose pathologies like "the EM in each
  # bootstrap iteration collapses to ran_eff_sd ~ 0", which would inflate the
  # prediction-MSE summary by removing the refit's ability to recover u_d*.
  boot_diag <- NULL
  if (!is.null(mse_estimated)) {
    qsum <- function(x) {
      x <- x[is.finite(x)]
      if (length(x) == 0L) return(rep(NA_real_, 5L))
      as.numeric(quantile(x, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE))
    }
    boot_diag <- list(
      ran_eff_sd_orig    = megb_fit$megb_model$ran_eff_sd,
      error_sd_orig      = megb_fit$megb_model$error_sd,
      boot_ran_eff_sd    = mse_estimated$boot_ran_eff_sd_boot,
      boot_error_sd      = mse_estimated$boot_error_sd,
      boot_ran_eff_sd_q  = qsum(mse_estimated$boot_ran_eff_sd_boot),
      boot_error_sd_q    = qsum(mse_estimated$boot_error_sd)
    )
    message(sprintf(
      "Bootstrap diagnostics: original ran_eff_sd = %.4f, error_sd = %.4f.\n  Bootstrap ran_eff_sd quartiles (min,Q1,Q2,Q3,max) = %s\n  Bootstrap error_sd   quartiles                  = %s\n  Frac bootstrap iters with ran_eff_sd < 1e-4 = %.2f",
      megb_fit$megb_model$ran_eff_sd %||% NA_real_,
      megb_fit$megb_model$error_sd   %||% NA_real_,
      paste(sprintf("%.4f", boot_diag$boot_ran_eff_sd_q), collapse = ", "),
      paste(sprintf("%.4f", boot_diag$boot_error_sd_q),   collapse = ", "),
      mean(boot_diag$boot_ran_eff_sd < 1e-4, na.rm = TRUE)
    ))
  }

  # ── 10. Assemble result ───────────────────────────────────────────────────────
  result <- list(
    ind            = ind,
    var            = var_df,
    CI             = ci_df,
    ind_bench      = ind_bench,
    var_bench      = var_bench_df,
    CI_bench       = ci_bench_df,
    yhat           = yhat,
    model          = model_ref,
    megb_model     = megb_fit$megb_model,
    gbm_engine     = gbm_engine,
    smp_data       = smp_data,
    out_call       = out_call,
    transformation = transformation,
    framework      = fwk,
    boot_diag      = boot_diag
  )

  class(result) <- c("megb", "xgb", "povmap")
  return(result)
}
