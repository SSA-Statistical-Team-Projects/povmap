#' This function evaluates the performance of a Gradient Boosting model
#' using the provided set of hyperparameters. It fits the model and 
#' computes the best evaluation score based on the engine type (XGBoost,
#' LightGBM, or CatBoost).
#'
#' @param params A list of hyperparameters for the Gradient Boosting model.
#' @param gbm_engine The engine used for the Gradient Boosting model. It can be 
#'        one of "xgboost", "lightgbm", or "catboost".
#' @param features_train A data frame or matrix containing the training features.
#' @param response_train A numeric vector containing the training response variable.
#' 
#' @return A list containing:
#' \item{Score}{The best evaluation score (e.g., RMSE) of the model.}
#' \item{params}{The hyperparameters used in the model.}
#' \item{best_iter}{The best iteration (round) of the model during training.}
#' 
#' @import dplyr
#' @import purrr

eval_params <- function(params, gbm_engine, features_train, response_train) {
  model <- train_gbmodel(
    model_type     = gbm_engine,
    features_train = features_train,
    response_train = response_train,
    params         = params
  )
  if (gbm_engine == "xgboost") {
    best_score <- min(model$eval_log[[grep("^test", names(model$eval_log), value = TRUE)[1]]])
  } else if (gbm_engine == "lightgbm") {
    best_score <- min(model$eval_log[, "valid"])
  } else if (gbm_engine == "catboost") {
    best_score <- min(model$eval_log[grep("^test", names(model$eval_log))[1]])
  }
  list(Score = as.numeric(best_score), params = unlist(params), best_iter = model$best_iter)
}

# -------------------------
# Funktion: Tuning Wrapper
# -------------------------
#' Tune Gradient Boosting Models with Hyperparameter Optimization
#'
#' This function tunes a Gradient Boosting model (XGBoost, LightGBM, or CatBoost)
#' using grid search, random search, or Bayesian optimization.
#'
#' @param gbm_engine A string indicating the model engine: "xgboost", "lightgbm", or "catboost".
#' @param Y A numeric vector of the response variable.
#' @param X A data frame of predictor variables.
#' @param smp_data A data frame with the training data.
#' @param param_space A list of parameter ranges for tuning.
#' @param tuning_strategy The strategy for tuning: "grid", "random", or "bayesian".
#' @param n_iter integer, number of hyperparameter configurations to be evaluated when using random search or Bayesian optimization
#' @param seed The random seed for reproducibility.
#'
#' @return A list containing the best parameters, the best score (RMSE), and a data frame with all tuning results.
#' \item{best_params}{The best hyperparameters found during the tuning.}
#' \item{best_score}{The best evaluation score (e.g., RMSE).}
#' \item{results_df}{A data frame containing the evaluation results for all hyperparameter combinations.}
#' \item{best_iter}{The best iteration of the model.}
#' \item{extra}{A list with additional information, such as the Bayesian Optimization object and stop status.}
#' @import dplyr
#' @import purrr
#' @import ParBayesianOptimization
#' @export


tune_gbmodel <- function(
    gbm_engine, 
    Y, 
    X,
    smp_data, 
    param_space,
    tuning_strategy = c("grid", "random", "bayesian"),
    n_iter = 20,
    seed = 42
) {
  set.seed(seed)
  tuning_strategy <- match.arg(tuning_strategy)
  
  if (gbm_engine %in% c("xgboost", "lightgbm")) {
    # identically to megb():
    mm_formula     <- terms(~ . - 1, data = smp_data[, names(X), drop = FALSE])
    features_train <- model.matrix(mm_formula, data = smp_data)
  } else if (gbm_engine == "catboost") {
    # CatBoost: native handling → no One-Hot-Encoding!
    features_train <- smp_data[, names(X), drop = FALSE]
  }
  response_train <- Y
  
  run_eval <- function(p) eval_params(p, gbm_engine, features_train, response_train)
  
  if (tuning_strategy == "grid") {
    param_grid <- expand.grid(param_space, stringsAsFactors = FALSE)
    results <- purrr::map(seq_len(nrow(param_grid)), ~run_eval(as.list(param_grid[.x, ])), .progress = TRUE)
  }
  
  if (tuning_strategy == "random") {
    sampler <- function(space) lapply(space, function(v) sample(v, 1))
    results <- purrr::map(seq_len(n_iter), ~run_eval(sampler(param_space)), .progress = TRUE)
  }
  
  if (tuning_strategy == "bayesian") {
    # --- encode discrete candidates in Index spaces ---
    enc_space  <- lapply(param_space, as.vector)
    idx_bounds <- setNames(lapply(enc_space, function(v) c(1L, length(v))), names(enc_space))
    
    # --- ParBayesianOptimization-Constraints ---
    n_inputs <- length(idx_bounds)
    if (n_iter <= n_inputs) {
      stop(sprintf(
        "Bayesian: n_iter (%d) must be > number of parameters (%d). Increase n_iter or reduce param_space.",
        n_iter, n_inputs
      ))
    }
    # initPoints > n_inputs and < n_iter; iters.n >= 1; iters.k = 1
    initPoints <- min(max(n_inputs + 1L, 5L), n_iter - 1L)
    iters.n    <- max(1L, n_iter - initPoints)
    iters.k    <- 1L
    
    # --- Logging all tries (for results_df) ---
    trial_log <- list()
    
    # --- Index -> real value (clamp & round) ---
    decode <- function(idx_args_named) {
      mapply(
        function(v, i) v[[ max(1L, min(length(v), round(as.numeric(i)))) ]],
        enc_space,
        idx_args_named[names(enc_space)],
        SIMPLIFY = FALSE
      )
    }
    
    # --- target function for BO: - maximise RMSE ---
    bayes_fun <- function(...) {
      idx_list <- list(...)                     # collect ingoing indices
      params <- decode(idx_list)                # translate indices to real values
      out <- run_eval(params)
      trial_log[[length(trial_log) + 1L]] <<- list(
        params    = out$params,
        Metric    = out$Score,     # loss metric (RMSE)
        best_iter = out$best_iter
      )
      list(Score = -out$Score)     # maximize internally
    }
    
    # --- BO execution ---
    opt <- bayesOpt(
      FUN        = bayes_fun,
      bounds     = idx_bounds,
      initPoints = initPoints,
      iters.n    = iters.n,
      iters.k    = iters.k,
      acq        = "ucb",
      plotProgress = F
    )
    results_df <- do.call(rbind, lapply(trial_log, function(z) {
      data.frame(t(as.data.frame(z$params)),
                 Metric    = z$Metric,
                 best_iter = z$best_iter,
                 check.names = FALSE)
    }))
    
    best_idx   <- which.min(results_df$Metric)
    best_params <- as.list(results_df[best_idx, setdiff(names(results_df), c("Metric","best_iter")), drop = FALSE])
    best_score  <- results_df$Metric[best_idx]
    best_iter   <- results_df$best_iter[best_idx]
    
    return(list(
      best_params = best_params,
      best_score  = best_score,
      results_df  = results_df,
      best_iter   = best_iter, 
      extra       = list(bo_object = opt, stop_status = opt$stopStatus)
    ))
  }
  
  # --- results as DataFrame ---
  df_results <- purrr::map_dfr(
    results,
    function(x) {
      data.frame(
        t(as.data.frame(x$params)),
        RMSE = x$Score,
        best_iter = x$best_iter,
        stringsAsFactors = FALSE
      )
    }
    
  )
  
  # --- find best result (always minimum RMSE) ---
  best_idx <- which.min(df_results$RMSE)
  best_params <- as.list(df_results[best_idx, setdiff(names(df_results), "RMSE")])
  best_score  <- df_results$RMSE[best_idx]
  best_iter   <- df_results$best_iter[best_idx]
  # --- Output ---
  list(
    best_params = best_params,
    best_score  = best_score,
    results_df  = df_results,
    best_iter   = best_iter
  )
}