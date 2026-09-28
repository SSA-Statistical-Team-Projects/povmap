# Internal: train a gradient boosting model for a given engine
#' @importFrom xgboost xgb.train xgb.cv xgb.DMatrix

train_gbmodel <- function(model_type, features_train, response_train,
                           params, cat_idx0 = integer(0),
                           weights = NULL, fit_full = TRUE, ...) {
  # weights: numeric vector of observation weights, length = nrow(features_train).
  # When NULL or all 1s, behaviour is unchanged from the historical default.
  #
  # fit_full (xgboost only): FALSE skips the full-data booster and returns the
  # cross-validation results (out-of-fold predictions, best_iter) plus the RNG
  # state at the point where the booster would have been trained. The one draw
  # xgb.train makes from R's RNG (its seed, when params carries none) is consumed
  # anyway, so the stream -- and every later fold assignment -- is unchanged.
  # em_gb_lmm uses this to train only the booster of the final EM iteration;
  # .xgb_full_fit() then reproduces that booster exactly. Requires out-of-fold
  # predictions; if they are unavailable the booster is trained as before.
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

      if (!fit_full && !is.null(oof)) {
        rng_state <- get(".Random.seed", envir = globalenv())
        if (!("seed" %in% names(xgb_params))) sample(.Machine$integer.max, size = 1)  # xgb.train's own draw
        return(list(engine         = "xgboost",
                    boosting       = NULL,
                    best_iter      = best_iter,
                    feature_names  = colnames(features_train),
                    prediction     = NULL,
                    oof_prediction = oof,
                    eval_log       = cv$evaluation_log,
                    train_pool_schema = NULL,
                    deferred       = list(rng_state = rng_state, params = xgb_params,
                                          nrounds = best_iter, weights = weights)))
      }

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

# Internal: train the full-data xgboost booster that train_gbmodel(fit_full = FALSE)
# skipped, exactly as train_gbmodel would have: same DMatrix construction, same
# params and rounds, and R's RNG set to the state at the skipped call so xgb.train
# draws the same seed. The caller's RNG state is restored afterwards.
.xgb_full_fit <- function(deferred, features_train, response_train) {
  dtrain <- xgboost::xgb.DMatrix(data  = data.matrix(features_train),
                                  label = response_train)
  if (!is.null(deferred$weights)) xgboost::setinfo(dtrain, "weight", as.numeric(deferred$weights))
  state_after <- get(".Random.seed", envir = globalenv())
  assign(".Random.seed", deferred$rng_state, envir = globalenv())
  on.exit(assign(".Random.seed", state_after, envir = globalenv()), add = TRUE)
  booster <- xgboost::xgb.train(data = dtrain, params = deferred$params,
                                nrounds = deferred$nrounds, verbose = 0)
  list(boosting = booster, prediction = stats::predict(booster, dtrain))
}
