# Internal: EM algorithm combining Gradient Boosting with a Linear Mixed Model
#' @importFrom lme4 lmer ranef VarCorr fixef
#' @importFrom xgboost xgboost xgb.cv
#' @importFrom stats logLik predict sigma

em_gb_lmm <- function(Y,
                      X,
                      formula_random_effects,
                      data,
                      gradient_params,
                      initial_random_effects,
                      max_iterations,
                      error_tolerance,
                      dom_name,
                      cov_names,
                      gbm_engine = "xgboost",
                      weights    = NULL,
                      ...) {
  # weights: numeric vector of observation weights aligned to rows of X / data.
  # Forwarded to train_gbmodel (xgb.DMatrix weight) AND to lme4::lmer (weights).
  # When NULL, behaviour is unchanged.
  #
  # Rescale survey weights WITHIN each domain so that, in every domain, the
  # weights sum to that domain's sample size (equivalently, average to 1). This
  # is the Pfeffermann et al. (1998) scaling used for survey-weighted multilevel
  # / pseudo-EBP estimation (Guadarrama et al.; Dang et al. 2026), and it
  # matches the rescaling the xgb pipeline applies (xgb.R: w * domain_n /
  # domain_sum_wts). It prevents domains with concentrated weights from
  # dominating the fit and keeps the weight scale unit-level for numerical
  # stability in lme4 / xgboost. `data[[dom_name]]` gives the domain of each row.
  if (!is.null(weights)) {
    dom   <- as.character(data[[dom_name]])
    dsum  <- ave(weights, dom, FUN = function(z) sum(z, na.rm = TRUE))
    dn    <- ave(weights, dom, FUN = length)
    keep  <- is.finite(dsum) & dsum > 0
    weights[keep]  <- weights[keep] * dn[keep] / dsum[keep]
    weights[!keep] <- 1
  }

  target            <- Y
  continue_condition <- TRUE
  iterations        <- 0
  dom_name_effects  <- 0
  old_log_lik       <- 0
  features_train    <- X

  # Initialise the EM with a LMM-only fit on Y so the first GB iteration sees
  # a domain-demeaned target rather than raw Y. Starting from zero random effects
  # lets GB absorb domain means in iteration 1, after which lme4 finds no residual
  # domain structure and the EM collapses to the trivial fixed point (ran_eff_sd=0).
  formula_lmm_init <- as.formula(paste0("target ~ 1 + ", formula_random_effects))
  # Attach weights as a column on `data` so lme4 can find them by name (works
  # around lme4's evaluation of `weights` in the model frame).
  if (!is.null(weights)) data$.megb_w <- as.numeric(weights)
  suppressMessages(
    lmefit_init <- lme4::lmer(
      formula_lmm_init, data = data, REML = FALSE,
      weights = if (!is.null(weights)) data$.megb_w else NULL
    )
  )
  adjusted_target <- target - (stats::predict(lmefit_init) - lme4::fixef(lmefit_init))

  while (continue_condition) {
    iterations <- iterations + 1

    response_train <- adjusted_target
    gbm_results    <- train_gbmodel(gbm_engine, features_train, response_train,
                                    params  = gradient_params,
                                    weights = weights)

    model        <- gbm_results$boosting
    unit_pred_smp <- gbm_results$prediction

    # Estimate the LMM (idiosyncratic variance sigma_e AND the shrunken random
    # effect) from OUT-OF-FOLD GB residuals, not in-sample residuals. A strong
    # learner overfits within-area in-sample, deflating sigma_e ~6x and inflating
    # the EBLUP shrinkage factor toward 1 (so even one-observation areas are
    # unshrunk and the RE over-absorbs covariate-predictable area structure). Using
    # held-out residuals restores the proper area-specific shrinkage gradient
    # (census-validated: gamma 0.99->~0.38 at n_d=1, ~0.85 at n_d>8) and removes the
    # over-absorption. The point prediction still uses the in-sample GB
    # (unit_pred_smp); only the variance/RE decomposition uses the held-out residual.
    # Applied every EM iteration so the corrected variance feeds the next iteration.
    # Disable with options(megb.oof_sigma = FALSE) to recover the legacy behaviour.
    .use_oof <- !isFALSE(getOption("megb.oof_sigma", TRUE)) &&
                !is.null(gbm_results$oof_prediction)
    tmp_res    <- if (.use_oof) Y - gbm_results$oof_prediction else Y - unit_pred_smp
    formula_lmm <- as.formula(paste0("tmp_res ~ 1 +", formula_random_effects))

    suppressMessages(
      lmefit <- lme4::lmer(
        formula_lmm, data = data, REML = FALSE,
        weights = if (!is.null(weights)) data$.megb_w else NULL
      )
    )

    new_log_lik <- as.numeric(stats::logLik(lmefit))

    continue_condition <- (
      abs((new_log_lik - old_log_lik[iterations]) / old_log_lik[iterations]) > error_tolerance &
        iterations < max_iterations
    )

    old_log_lik     <- c(old_log_lik, new_log_lik)
    dom_name_effects <- stats::predict(lmefit)
    adjusted_target  <- target - stats::predict(lmefit) + lme4::fixef(lmefit)
  }

  residuals  <- target - stats::predict(lmefit) - unit_pred_smp
  error_sd   <- stats::sigma(lmefit)

  importance_matrix <- get_feature_importance(
    engine         = gbm_results$engine,
    model          = gbm_results$boosting,
    features_train = features_train,
    train_pool     = gbm_results$train_pool_schema
  )

  list(
    call                 = match.call(),
    boosting             = model,
    effect_model         = lmefit,
    dom_name_effects     = lme4::ranef(lmefit),
    ran_eff_sd           = as.data.frame(lme4::VarCorr(lmefit))$sdcor[1],
    error_sd             = error_sd,
    variance_covariance  = lme4::VarCorr(lmefit),
    log_lik              = old_log_lik,
    iterations_used      = iterations,
    residuals            = residuals,
    dom_name             = dom_name,
    initial_random_effects = initial_random_effects,
    importance_matrix    = importance_matrix,
    eval_log             = gbm_results$eval_log,
    gradient_params      = gradient_params
  )
}
