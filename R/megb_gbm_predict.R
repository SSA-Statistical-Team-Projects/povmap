# Internal: generate unit-level predictions from a fitted MEGB model
#' @importFrom lme4 fixef
#' @importFrom stats predict

gbm_predict <- function(model, smp_data, pop_data, Y, dom_name, gbm_engine, cov_names) {

  effect_model <- model$effect_model
  re_smp <- stats::predict(effect_model, smp_data, allow.new.levels = TRUE)
  re_pop <- stats::predict(effect_model, pop_data, allow.new.levels = TRUE)
  fe     <- lme4::fixef(effect_model)

  switch(
    gbm_engine,

    "xgboost" = {
      x_smp <- data.matrix(smp_data[, cov_names, drop = FALSE])
      x_pop <- data.matrix(pop_data[, cov_names, drop = FALSE])
      gb_smp <- stats::predict(model$boosting, x_smp)
      gb_pop <- stats::predict(model$boosting, x_pop)
      unit_pred_smp <- gb_smp + re_smp - fe
      res           <- as.numeric(Y) - unit_pred_smp
      unit_preds    <- gb_pop + re_pop - fe
      list(unit_pred_smp = unit_pred_smp, res = res,
           unit_pred_pop = data.frame(unit_preds = unit_preds,
                                      dom_name   = pop_data[[dom_name]],
                                      row.names  = NULL))
    },

    "lightgbm" = {
      x_smp <- data.matrix(smp_data[, cov_names, drop = FALSE])
      x_pop <- data.matrix(pop_data[, cov_names, drop = FALSE])
      gb_smp <- stats::predict(model$boosting, x_smp)
      gb_pop <- stats::predict(model$boosting, x_pop)
      unit_pred_smp <- gb_smp + re_smp - fe
      res           <- as.numeric(Y) - unit_pred_smp
      unit_preds    <- gb_pop + re_pop - fe
      list(unit_pred_smp = unit_pred_smp, res = res,
           unit_pred_pop = data.frame(unit_preds = unit_preds,
                                      dom_name   = pop_data[[dom_name]],
                                      row.names  = NULL))
    },

    "catboost" = {
      df_smp   <- smp_data[, cov_names, drop = FALSE]
      df_pop   <- pop_data[, cov_names, drop = FALSE]
      pool_smp <- catboost::catboost.load_pool(data = df_smp)
      pool_pop <- catboost::catboost.load_pool(data = df_pop)
      gb_smp   <- as.numeric(catboost::catboost.predict(model$boosting, pool_smp))
      gb_pop   <- as.numeric(catboost::catboost.predict(model$boosting, pool_pop))
      unit_pred_smp <- gb_smp + re_smp - fe
      res           <- as.numeric(Y) - unit_pred_smp
      unit_preds    <- gb_pop + re_pop - fe
      list(unit_pred_smp = unit_pred_smp, res = res,
           unit_pred_pop = data.frame(unit_preds = unit_preds,
                                      dom_name   = pop_data[[dom_name]],
                                      row.names  = NULL))
    }
  )
}
