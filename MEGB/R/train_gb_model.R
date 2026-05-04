#' train GB model - MEGB helper function
#' @param model_type character, indicating which gbm_engine was used
#' @param features_train covariates used for training
#' @param response_train response var used for training
#' @param params parameters to use for models
#' @param cat_idx0 index for catboost, default is 0
#' @param ... additional arguments
#' @import xgboost
#' @import lightgbm
#' @import catboost
#' @return list containing boosting information

train_gbmodel <- function(model_type, 
                          features_train, 
                          response_train, 
                          params, 
                          cat_idx0 = integer(0), 
                          ...) {
  
  model <- switch(
    model_type,
    
    # --- XGBoost ---
    "xgboost" = {
      nrounds    <- params$nrounds %||% 1000
      xgb_params <- params[setdiff(names(params), "nrounds")]
      
      cv <- xgb.cv(
        params  = xgb_params,
        data    = data.matrix(features_train),
        label   = response_train,
        nrounds = nrounds,
        nfold   = 5,
        early_stopping_rounds = 10,
        verbose = 0,
        prediction = FALSE,
        ...
      )
      
      best_iter <- cv$best_iteration
      
      booster <- xgboost(
        data    = data.matrix(features_train),
        label   = response_train,
        params  = xgb_params,
        nrounds = best_iter,
        verbose = 0
      )
      
      list(
        engine        = "xgboost",
        boosting      = booster,
        best_iter     = best_iter,
        feature_names = colnames(features_train),
        prediction    = predict(booster, data.matrix(features_train)),
        eval_log      = cv$evaluation_log   
      )
    },
    
    # --- LightGBM ---
    "lightgbm" = {
      nrounds    <- params$nrounds %||% 1000
      dtrain <- lgb.Dataset(data.matrix(features_train), label = response_train)
      
      cv <- lgb.cv(
        params  = params,
        data    = dtrain,
        nrounds = nrounds,
        nfold   = 5,
        early_stopping_rounds = 10,
        verbose = -1,
        eval_train_metric = TRUE,
        ...
      )
      
      best_iter <- cv$best_iter
      
      booster <- lgb.train(
        params  = params,
        data    = dtrain,
        nrounds = best_iter,
        verbose = -1
      )
      
      # For predictions with LightGBM Matrix, not Dataset:
      pred_train <- predict(booster, data.matrix(features_train))
      
      # Eval Log 
      train_rmse <- unlist(cv$record_evals$train$rmse$eval)
      valid_rmse <- unlist(cv$record_evals$valid$rmse$eval)
      eval_log <- data.frame(
        iter  = seq_along(train_rmse),
        train = as.numeric(train_rmse),
        valid = as.numeric(valid_rmse)
      )
      
      
      list(
        engine        = "lightgbm",
        boosting      = booster,
        best_iter     = best_iter,
        feature_names = colnames(features_train),
        prediction    = pred_train,
        eval_log      = eval_log 
      )
    },
    
    # --- CatBoost ---
    "catboost" = {
      train_pool <- catboost.load_pool(
        data  = features_train,
        label = response_train
      )
      
      cv <- catboost.cv(
        pool   = train_pool,
        params = params,
        early_stopping_rounds = 10,
        ...
      )
      
      # select first numeric metric column (z. B. test-RMSE)
      numeric_cols <- vapply(cv, is.numeric, logical(1))
      metric_col   <- names(cv)[numeric_cols][1]
      best_iter    <- which.min(cv[[metric_col]])
      
      params$iterations <- best_iter
      
      booster <- catboost.train(
        learn_pool = train_pool,
        params     = params
      )
      
      pred_train <- as.numeric(catboost.predict(booster, train_pool))
      
      # Eval-Log: catboost.cv already returns a Data Frame of the metrics
      eval_log <- cv
      if (!("iter" %in% names(eval_log))) {
        eval_log$iter <- seq_len(nrow(eval_log))
      }
      
      list(
        engine            = "catboost",
        boosting          = booster,
        best_iter         = best_iter,
        feature_names     = colnames(features_train),
        train_pool_schema = train_pool,
        prediction        = pred_train,
        eval_log          = eval_log       
      )
    }
  )
  
  model
}

