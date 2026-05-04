#' MEGB helper function
#' @param model model object
#' @param smp_data dataframe, sample data
#' @param pop_data dataframe, population data
#' @param Y target variable
#' @param dom_name factor, domain name
#' @param gbm_engine character, name of gbm engine to use
#' @param cov_names names of covariates
#' @import lme4
#' @import stats
#' @import xgboost
#' @import lightgbm
#' @import catboost
#' @return list of unit predictions

gbm_predict <- function(model,
                        smp_data,
                        pop_data,
                        Y,
                        dom_name,
                        gbm_engine,
                        cov_names
                        ) {

  effect_model <- model$effect_model
  
  # RE/Fixef-part preparation
  re_smp <- predict(effect_model, smp_data, allow.new.levels = TRUE)
  re_pop <- predict(effect_model, pop_data, allow.new.levels = TRUE)
  fe     <- fixef(effect_model)
  
  switch(
    gbm_engine,
    
    # --- XGBoost ---
    "xgboost" = {
      x_smp <- data.matrix(smp_data[, cov_names, drop = FALSE])
      x_pop <- data.matrix(pop_data[, cov_names, drop = FALSE])
      
      gb_smp <- predict(model$boosting, x_smp)
      gb_pop <- predict(model$boosting, x_pop)
      
      unit_pred_smp <- gb_smp + re_smp - fe
      res           <- as.numeric(Y) - unit_pred_smp
      unit_preds    <- gb_pop + re_pop - fe
      
      list(unit_pred_smp = unit_pred_smp,
           res            = res,
           unit_pred_pop  = data.frame(unit_preds = unit_preds,
                                       dom_name   = pop_data[[dom_name]],
                                       row.names  = NULL))
    },
    
    # --- LightGBM ---
    "lightgbm" = {
      x_smp <- data.matrix(smp_data[, cov_names, drop = FALSE])
      x_pop <- data.matrix(pop_data[, cov_names, drop = FALSE])
      
      gb_smp <- stats::predict(model$boosting, x_smp)
      gb_pop <- stats::predict(model$boosting, x_pop)
      
      unit_pred_smp <- gb_smp + re_smp - fe
      res           <- as.numeric(Y) - unit_pred_smp
      unit_preds    <- gb_pop + re_pop - fe
      
      list(unit_pred_smp = unit_pred_smp,
           res            = res,
           unit_pred_pop  = data.frame(unit_preds = unit_preds,
                                       dom_name   = pop_data[[dom_name]],
                                       row.names  = NULL))
    },
    
    # --- CatBoost ---
    "catboost" = {
      # Keep data.frame to preserve factor columns; DO NOT coerce to matrix
      df_smp <- smp_data[, cov_names, drop = FALSE]
      df_pop <- pop_data[, cov_names, drop = FALSE]
      
      pool_smp <- catboost.load_pool(data         = df_smp)
      pool_pop <- catboost.load_pool(data         = df_pop)
      
      gb_smp <- as.numeric(catboost.predict(model$boosting, pool_smp))
      gb_pop <- as.numeric(catboost.predict(model$boosting, pool_pop))
      
      gb_smp <- as.numeric(catboost.predict(model$boosting, pool_smp))
      gb_pop <- as.numeric(catboost.predict(model$boosting, pool_pop))
      
      unit_pred_smp <- gb_smp + re_smp - fe
      res           <- as.numeric(Y) - unit_pred_smp
      unit_preds    <- gb_pop + re_pop - fe
      
      list(unit_pred_smp = unit_pred_smp,
           res            = res,
           unit_pred_pop  = data.frame(unit_preds = unit_preds,
                                       dom_name   = pop_data[[dom_name]],
                                       row.names  = NULL))
    }
  )

}
