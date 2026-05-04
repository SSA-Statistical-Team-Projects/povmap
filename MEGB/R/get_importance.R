#' get feature importance - MEGB helper function
#' @param engine character, used engine
#' @param model list, megb model
#' @param features_train covariates used for training, default is NULL
#' @param train_pool train pool for catboost feature importance estimation, default is NULL
#' @import xgboost
#' @import lightgbm
#' @import catboost

get_feature_importance <- function(engine, model, features_train = NULL, train_pool = NULL) {
  switch(engine,
         
         "xgboost" = {
           xgb.importance(model = model)
         },
         
         "lightgbm" = {
           lgb.importance(model, percentage = TRUE)
         },
         
         "catboost" = {
           imp <- catboost.get_feature_importance(
             model = model,
             pool  = train_pool,
             type  = "FeatureImportance"
           )
           data.frame(
             Feature   = colnames(features_train),
             Importance = imp
           )
         },
         
         stop("Engine not supported.")
  )
}
