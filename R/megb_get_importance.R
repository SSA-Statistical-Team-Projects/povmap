# Internal: compute feature importance for a fitted GB model
#' @importFrom xgboost xgb.importance

get_feature_importance <- function(engine, model, features_train = NULL,
                                    train_pool = NULL) {
  switch(
    engine,
    "xgboost"  = xgboost::xgb.importance(model = model),
    "lightgbm" = {
      if (!requireNamespace("lightgbm", quietly = TRUE))
        return(data.frame(Feature = character(0), Importance = numeric(0)))
      lightgbm::lgb.importance(model, percentage = TRUE)
    },
    "catboost" = {
      if (!requireNamespace("catboost", quietly = TRUE))
        return(data.frame(Feature = character(0), Importance = numeric(0)))
      imp <- catboost::catboost.get_feature_importance(
        model = model, pool = train_pool, type = "FeatureImportance"
      )
      data.frame(Feature = colnames(features_train), Importance = imp)
    },
    stop("Engine not supported.")
  )
}
