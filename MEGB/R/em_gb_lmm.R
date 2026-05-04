
#' EM-algorithm for Mixed Effect Gradient Boosting
#'
#' This function combines Gradient Boosting with a Linear Mixed
#' Model (LMM) by using the Expectation-Maximization (EM) algorithm.
#' @param Y numeric, target variable
#' @param X data.frame, Covariates
#' @param formula_random_effects formula, formula for Random Effects
#' @param data data.frame
#' @param gradient_params list, optional list of Gradient Boosting parameters
#' @param initial_random_effects numeric, initial value for domain effects
#' @param max_iterations integer, max number of EM-iterations
#' @param error_tolerance numeric, error tolerance for stop-criterion
#' @param dom_name factor, name of domain variable
#' @param cov_names name of covariates
#' @param gbm_engine character, choice of Gradient Boosting implementation ("xgboost", "lightgbm", "catboost")
#' @param ... additional arguments
#' @import magrittr
#' @import pbapply
#' @import xgboost
#' @import catboost
#' @import lightgbm
#' @importFrom lme4 lmer ranef VarCorr fixef
#' @return list with results

em_gb_lmm <- function(Y,                          # target variable
                      X,                          # covariates
                      formula_random_effects,     # formula for Random Effects
                      data,                       # data
                      gradient_params,            # Gradient Boosting parameters
                      initial_random_effects,     # initial values for domain-effects
                      max_iterations,             # max EM-iterations
                      error_tolerance,            # error tolerance for stop-criterion
                      dom_name,                   # name of domain variable
                      cov_names,                  # names of covariates,
                      gbm_engine = "xgboost",     # choice of Gradient Boosting implementation ("xgboost", "lightgbm", "catboost")
                      ...) {                      # additional arguments

  # Initialisation
  target <- Y                                      # target variable
  continue_condition <- TRUE                      # looping condition
  iterations <- 0                                 # counter for iterations
  dom_name_effects <- 0                           # initial values for domain-effects
  adjusted_target <- target - initial_random_effects # adjusted target variable
  old_log_lik <- 0                                # starting value for Log-Likelihood
  rmse <- NULL                                    # Root Mean Square Error (optional)
  features_train <- X
  
  # EM-algorithm loop
  while (continue_condition) {
    iterations <- iterations + 1

    # Gradient Boosting -----------------------------------------------------

    response_train <- adjusted_target
    
    gbm_results <- train_gbmodel(gbm_engine,
                                   features_train,
                                   response_train,
                                   params = gradient_params)
    
    model <- gbm_results$boosting
    unit_pred_smp <- gbm_results$prediction
    best_iter <- gbm_results$best_iter
    
    # if (gbm_engine %in% c("xgboost", "lightgbm")) {
    #   gradient_params$nrounds   <- best_iter
    #   gradient_params$iterations <- NULL   # remove for safety
    # } else if (gbm_engine == "catboost") {
    #   gradient_params$iterations <- best_iter
    #   gradient_params$nrounds    <- NULL   # remove for safety
    # }
    
    
    # calculate residuals for LMM
    tmp_res <- Y - unit_pred_smp

    # Linear Mixed Model (LMM) --------------------------------------
    # formula for LMM with Random Effects
    formula_lmm <- as.formula(paste0("tmp_res ~ 1 +", formula_random_effects))

    suppressMessages(
    # fit LMM to data
    lmefit <- lme4::lmer(formula_lmm, data = data, REML = FALSE)
    )
    
    # save Log-Likelihood of LMM
    new_log_lik <- as.numeric(stats::logLik(lmefit))

    # EM algorithm --------------------------------------------------------
    # stop-criterion: Changes in Log-Likelihood or max Iterations
    continue_condition <- (
      abs((new_log_lik - old_log_lik[iterations]) / old_log_lik[iterations]) > error_tolerance &
        iterations < max_iterations
    )

    # renew Log-Likelihood
    old_log_lik <- c(old_log_lik, new_log_lik)

    # new domain effects from the LMM
    dom_name_effects <- predict(lmefit) 

    # adjust target variable: old target variable minus LMM-Predictions
    adjusted_target <- target - predict(lmefit) + fixef(lmefit)
  }

  
  # calculate final residuals
  residuals <- target - predict(lmefit) - unit_pred_smp
  error_sd <- stats::sigma(lmefit)
  # calculate Feature Importance
  importance_matrix <- get_feature_importance(
    engine        = gbm_results$engine,
    model         = gbm_results$boosting,
    features_train = features_train,
    train_pool    = gbm_results$train_pool_schema
  )
  
  # save results
  result <- list(
    call = call,
    boosting = model,                                # Gradient Boosting model
    effect_model = lmefit,                           # LMM model
    dom_name_effects = lme4::ranef(lmefit),          # Random Effects
    ran_eff_sd = as.data.frame(lme4::VarCorr(lmefit))$sdcor[1], # standard deviation from Random Effects
    error_sd = error_sd,                             # standard deviation of errors
    variance_covariance = lme4::VarCorr(lmefit),     # variance-covariance-matrix
    log_lik = old_log_lik,                           # Log-Likelihood progression
    iterations_used = iterations,                    # number of iterations
    residuals = residuals,                           # residuals
    dom_name = dom_name,                             # Name der Domainvariable
    initial_random_effects = initial_random_effects, # initial value of random effects
    importance_matrix = importance_matrix,           # importance matrix
    eval_log = gbm_results$eval_log,                 # evaluation log for loss visualisation,
    gradient_params = gradient_params                # updated GB parameters
  )

  # return results
  result
}

