# Internal: train a gradient boosting model for a given engine
#' @importFrom xgboost xgb.train xgb.cv xgb.DMatrix

train_gbmodel <- function(model_type, features_train, response_train,
                           params, cat_idx0 = integer(0),
                           weights = NULL, ...) {
  # weights: numeric vector of observation weights, length = nrow(features_train).
  # When NULL or all 1s, behaviour is unchanged from the historical default.
  switch(
    model_type,

    "xgboost" = {
      nrounds    <- if (!is.null(params$nrounds)) params$nrounds
                   else if (!is.null(params$nround)) params$nround
                   else 1000
      xgb_params <- params[setdiff(names(params), c("nrounds", "nround"))]
      # Strip NULL / zero-length entries before they reach xgboost. Recent
      # xgboost versions serialise an empty param to '{}' and reject it with
      # "Invalid Parameter format for reg_alpha expect float but value='{}'",
      # which surfaces here whenever a tune slot (e.g. all_tunes[[k]]$alpha)
      # was saved as numeric(0) or NULL.
      xgb_params <- xgb_params[lengths(xgb_params) > 0L]

      dtrain <- xgboost::xgb.DMatrix(data  = data.matrix(features_train),
                                      label = response_train)
      if (!is.null(weights)) xgboost::setinfo(dtrain, "weight", as.numeric(weights))

      cv <- xgboost::xgb.cv(
        params                = xgb_params,
        data                  = dtrain,
        nrounds               = nrounds,
        nfold                 = 5,
        early_stopping_rounds = 10,
        verbose               = 0,
        prediction            = TRUE,   # keep out-of-fold preds for honest sigma_e
        ...
      )
      best_iter <- cv$best_iteration
      if (is.null(best_iter) || length(best_iter) == 0 || best_iter <= 0)
        best_iter <- nrounds
      # Out-of-fold (held-out) predictions, used by em_gb_lmm to estimate the
      # idiosyncratic error variance from held-out rather than in-sample residuals.
      # A strong learner overfits within-area in-sample, which deflates sigma_e and
      # inflates the EBLUP shrinkage toward 1 (the RE then over-absorbs covariate-
      # predictable area structure). xgboost >= 2.x returns these in $cv_predict;
      # older versions in $pred.
      oof <- cv$cv_predict; if (is.null(oof)) oof <- cv$pred
      if (is.list(oof)) oof <- oof[[1]]
      oof <- suppressWarnings(as.numeric(oof))
      if (length(oof) != nrow(features_train)) oof <- NULL  # guard: fall back to in-sample

      booster <- xgboost::xgb.train(
        data    = dtrain,
        params  = xgb_params,
        nrounds = best_iter,
        verbose = 0
      )

      list(engine        = "xgboost",
           boosting      = booster,
           best_iter     = best_iter,
           feature_names = colnames(features_train),
           prediction    = stats::predict(booster, dtrain),
           oof_prediction = oof,
           eval_log      = cv$evaluation_log,
           train_pool_schema = NULL)
    },

    "lightgbm" = {
      if (!requireNamespace("lightgbm", quietly = TRUE))
        stop("Package 'lightgbm' is required for gbm_engine = 'lightgbm'.")
      nrounds <- if (!is.null(params$nrounds)) params$nrounds else 1000
      dtrain  <- lightgbm::lgb.Dataset(data.matrix(features_train), label = response_train,
                                       weight = if (!is.null(weights)) as.numeric(weights) else NULL)

      cv <- lightgbm::lgb.cv(
        params                = params,
        data                  = dtrain,
        nrounds               = nrounds,
        nfold                 = 5,
        early_stopping_rounds = 10,
        verbose               = -1,
        eval_train_metric     = TRUE,
        ...
      )
      best_iter <- cv$best_iter

      booster <- lightgbm::lgb.train(
        params  = params,
        data    = dtrain,
        nrounds = best_iter,
        verbose = -1
      )
      pred_train  <- stats::predict(booster, data.matrix(features_train))
      train_rmse  <- unlist(cv$record_evals$train$rmse$eval)
      valid_rmse  <- unlist(cv$record_evals$valid$rmse$eval)
      eval_log    <- data.frame(iter  = seq_along(train_rmse),
                                train = as.numeric(train_rmse),
                                valid = as.numeric(valid_rmse))

      list(engine        = "lightgbm",
           boosting      = booster,
           best_iter     = best_iter,
           feature_names = colnames(features_train),
           prediction    = pred_train,
           eval_log      = eval_log,
           train_pool_schema = NULL)
    },

    "catboost" = {
      if (!requireNamespace("catboost", quietly = TRUE))
        stop("Package 'catboost' is required for gbm_engine = 'catboost'.")
      train_pool <- catboost::catboost.load_pool(data = features_train, label = response_train,
                                                 weight = if (!is.null(weights)) as.numeric(weights) else NULL)

      cv <- catboost::catboost.cv(pool = train_pool, params = params,
                                   early_stopping_rounds = 10, ...)

      numeric_cols <- vapply(cv, is.numeric, logical(1))
      metric_col   <- names(cv)[numeric_cols][1]
      best_iter    <- which.min(cv[[metric_col]])
      params$iterations <- best_iter

      booster    <- catboost::catboost.train(learn_pool = train_pool, params = params)
      pred_train <- as.numeric(catboost::catboost.predict(booster, train_pool))
      eval_log   <- cv
      if (!("iter" %in% names(eval_log))) eval_log$iter <- seq_len(nrow(eval_log))

      list(engine            = "catboost",
           boosting          = booster,
           best_iter         = best_iter,
           feature_names     = colnames(features_train),
           prediction        = pred_train,
           eval_log          = eval_log,
           train_pool_schema = train_pool)
    },

    stop("gbm_engine must be one of 'xgboost', 'lightgbm', or 'catboost'.")
  )
}
