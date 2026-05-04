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
#' @param mse logical. If \code{TRUE} (default), estimates MSE via the MEGB
#'   parametric bootstrap. Confidence intervals use a normal approximation with
#'   a delta-method correction when a transformation is applied.
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
                 cpus             = NULL,
                 ...) {

  out_call <- match.call()

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
  add_benchmark_megb <- function(x, benchmark_level, fwk, fixed,
                                 benchmark, benchmark_type) {
    point_estim       <- list(ind = data.frame(Mean = x))
    if (is.null(benchmark_level)) {
      point_estim$ind <- benchmark_ebp_national(
        point_estim    = point_estim,
        framework      = fwk,
        fixed          = fixed,
        benchmark      = benchmark,
        benchmark_type = benchmark_type)
    } else {
      point_estim$ind <- benchmark_ebp_level(
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

  # ── 6b. Benchmarking of point estimates ───────────────────────────────────────
  ind_bench <- NULL
  if (!is.null(benchmark)) {
    bench_vals <- add_benchmark_megb(
      x               = ind$Mean,
      benchmark_level = benchmark_level,
      fwk             = fwk,
      fixed           = fixed,
      benchmark       = benchmark,
      benchmark_type  = benchmark_type
    )
    ind_bench <- data.frame(
      Domain = ind$Domain,
      Mean   = bench_vals,
      stringsAsFactors = FALSE
    )
  }

  # ── 6c. Build benchmark_fn closure for the bootstrap ─────────────────────────
  if (!is.null(benchmark) && mse) {
    # Use a lean captured environment: explicit values only, parent = package ns.
    # This prevents serialising the full megb() scope to each parallel worker.
    benchmark_fn <- function(est_orig, boot_smp_data) {
      bm_val <- if (is.numeric(.bm) || is.data.frame(.bm)) {
        .bm
      } else if (is.null(.bm_level)) {
        c(Mean = mean(.corr_bt(boot_smp_data$y_star)))
      } else {
        y_bt <- .corr_bt(boot_smp_data$y_star)
        grp  <- boot_smp_data[[.bm_level]]
        lev_means <- tapply(y_bt, grp, mean)
        df <- data.frame(as.character(names(lev_means)),
                         as.numeric(lev_means),
                         stringsAsFactors = FALSE)
        names(df) <- c(.bm_level, "Mean")
        df
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
        .fwk      = fwk,
        .fixed    = fixed,
        .corr_bt  = corrected_bt,
        .add_bm   = add_benchmark_megb
      ),
      parent = getNamespace("povmap")
    )
  } else {
    benchmark_fn <- NULL
  }

  # ── 6d. Parametric bootstrap MSE (now that corrected_bt/benchmark_fn are ready)
  # pop_data_proc is one-hot encoded and lacks any benchmark_level column.
  # Re-attach it from fwk$pop_data so bootstrap samples carry it.
  boot_pop_data <- megb_fit$pop_data_proc
  if (!is.null(benchmark_level) && mse &&
      !benchmark_level %in% colnames(boot_pop_data)) {
    bm_col <- fwk$pop_data[[benchmark_level]]
    if (!is.null(bm_col) && length(bm_col) == nrow(boot_pop_data))
      boot_pop_data[[benchmark_level]] <- bm_col
    else
      warning("benchmark_level '", benchmark_level,
              "' not found in pop_data or row count mismatch — ",
              "benchmarked bootstrap MSE will be skipped.")
  }

  mse_estimated <- NULL
  if (mse) {
    message("Bootstrap with ", B, " iterations has started")
    mse_estimated <- mse_megb(
      Y                      = megb_fit$inp_smp_data$target_var,
      X                      = megb_fit$X_proc,
      dom_name               = domains,
      smp_data               = megb_fit$inp_smp_data$smp_data,
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
      unit_preds             = megb_fit$unit_preds_all,
      corrected_bt           = corrected_bt,
      benchmark_fn           = benchmark_fn
    )
  }

  # ── 7. MSE, delta-method variance correction, and CI ─────────────────────────
  var_df       <- NULL
  ci_df        <- NULL
  var_bench_df <- NULL
  ci_bench_df  <- NULL

  if (mse && !is.null(mse_estimated)) {

    z_val <- qnorm(1 - (1 - conf_level) / 2)

    # ── 7a. Unbenchmarked MSE (delta-method from transformed-space bootstrap) ──
    mse_raw        <- mse_estimated$MSE_estimates
    colnames(mse_raw)[1] <- "Domain"
    mse_raw$Domain <- as.character(mse_raw$Domain)
    colnames(mse_raw)[2] <- "MSE_t"

    mean_t        <- as.data.frame(megb_fit$Indicators)  # domain means in transformed space
    colnames(mean_t)[colnames(mean_t) == "dom_name"] <- "Domain"
    colnames(mean_t)[colnames(mean_t) == "Mean"]     <- "Mean_t"
    mean_t$Domain <- as.character(mean_t$Domain)

    ind$Domain <- as.character(ind$Domain)
    merged     <- merge(ind,    mse_raw, by = "Domain")
    merged     <- merge(merged, mean_t,  by = "Domain")

    # Delta-method: Var(g^{-1}(theta_hat_t)) ≈ MSE_t * [d g^{-1}/dx]^2
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

    if (transformation %in% c("arcsin", "logistic")) {
      ci_df$Lower <- pmax(ci_df$Lower, 0); ci_df$Upper <- pmin(ci_df$Upper, 1)
    } else if (transformation %in% c("poisson", "log", "log.shift")) {
      ci_df$Lower <- pmax(ci_df$Lower, 0)
    }

    ind <- data.frame(Domain = merged$Domain, Mean = merged$Mean,
                      stringsAsFactors = FALSE)

    # ── 7b. Benchmarked MSE (bootstrap computed in original space, no delta) ──
    if (!is.null(benchmark) && !is.null(mse_estimated$MSE_bench_estimates)) {

      mse_bench_raw        <- mse_estimated$MSE_bench_estimates
      colnames(mse_bench_raw)[1] <- "Domain"
      mse_bench_raw$Domain <- as.character(mse_bench_raw$Domain)
      colnames(mse_bench_raw)[2] <- "MSE_bench"

      ind_bench$Domain <- as.character(ind_bench$Domain)
      merged_b <- merge(ind_bench, mse_bench_raw, by = "Domain")

      var_bench_df <- data.frame(Domain = merged_b$Domain, Mean = merged_b$MSE_bench,
                                 stringsAsFactors = FALSE)

      ci_bench_df <- data.frame(
        Domain = merged_b$Domain,
        Lower  = merged_b$Mean - z_val * sqrt(pmax(merged_b$MSE_bench, 0)),
        Upper  = merged_b$Mean + z_val * sqrt(pmax(merged_b$MSE_bench, 0)),
        stringsAsFactors = FALSE
      )

      if (transformation %in% c("arcsin", "logistic")) {
        ci_bench_df$Lower <- pmax(ci_bench_df$Lower, 0)
        ci_bench_df$Upper <- pmin(ci_bench_df$Upper, 1)
      } else if (transformation %in% c("poisson", "log", "log.shift")) {
        ci_bench_df$Lower <- pmax(ci_bench_df$Lower, 0)
      }

      ind_bench <- data.frame(Domain = merged_b$Domain, Mean = merged_b$Mean,
                               stringsAsFactors = FALSE)
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
    framework      = fwk
  )

  class(result) <- c("megb", "xgb", "povmap")
  return(result)
}
