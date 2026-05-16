# Internal input validation for megb_em
#' @importFrom checkmate assert_numeric assert_data_frame assert_names
#' @importFrom checkmate assert_factor assert_character assert_list
#' @importFrom checkmate assert_flag assert_integerish assert_choice

input_checks_megb <- function(Y, X, dom_name, smp_data, pop_data,
                               gradient_params, na.rm, seed,
                               bootstrap_cores, mse, B, gbm_engine) {

  checkmate::assert_numeric(Y, any.missing = FALSE, min.len = 1,
                             finite = TRUE, null.ok = FALSE, len = nrow(X))
  # Allow NAs in X: xgboost/lightgbm/catboost all handle missing predictors
  # natively, and na.rm = TRUE filters incomplete rows downstream. Matches
  # xgb()'s permissive behaviour.
  checkmate::assert_data_frame(X, col.names = "named", min.rows = 1,
                                min.cols = 1, any.missing = TRUE)
  checkmate::assert_character(names(X), any.missing = FALSE, min.len = 1)
  checkmate::assert_data_frame(pop_data, min.rows = 1, any.missing = TRUE)
  checkmate::assert_data_frame(smp_data, min.rows = 1, any.missing = TRUE)
  checkmate::assert_names(intersect(colnames(pop_data), colnames(smp_data)),
                           must.include = c(dom_name))

  for (col_name in names(smp_data)) {
    if (is.character(smp_data[[col_name]])) {
      smp_data[[col_name]] <- as.factor(smp_data[[col_name]])
      warning("Variable ", col_name, " in smp_data coerced from character to factor.")
    }
  }
  checkmate::assert_factor(smp_data[[dom_name]], any.missing = FALSE)

  for (col_name in names(pop_data)) {
    if (is.character(pop_data[[col_name]])) {
      pop_data[[col_name]] <- as.factor(pop_data[[col_name]])
      warning("Variable ", col_name, " in pop_data coerced from character to factor.")
    }
  }
  checkmate::assert_factor(pop_data[[dom_name]], any.missing = FALSE)

  Y_name    <- tail(strsplit(deparse(substitute(Y)), "\\$")[[1]], 1)
  smp_vars  <- colnames(smp_data)[colnames(smp_data) != Y_name]
  missing_vars <- setdiff(smp_vars, colnames(pop_data))
  if (length(missing_vars) > 0)
    stop("Variables from smp_data missing in pop_data: ",
         paste(missing_vars, collapse = ", "))

  if (!all(unique(smp_data[[dom_name]]) %in% unique(pop_data[[dom_name]])))
    stop("There are domains in smp_data not present in pop_data.")

  checkmate::assert_choice(gbm_engine, choices = c("xgboost", "lightgbm", "catboost"))

  if (gbm_engine == "catboost" && !requireNamespace("catboost", quietly = TRUE))
    stop("Package 'catboost' is required for gbm_engine = 'catboost'.")

  if (is.null(gradient_params)) {
    message("gradient_params is NULL — default parameters used for ", gbm_engine, ".")
    gradient_params <- switch(gbm_engine,
      "xgboost"  = list(eta = 0.1, max_depth = 3, nrounds = 100, subsample = 1),
      "lightgbm" = list(objective = "regression", metric = "rmse",
                        learning_rate = 0.1, num_leaves = 63,
                        n_estimators = 100, bagging_fraction = 1, nrounds = 100),
      "catboost" = list(loss_function = "RMSE", learning_rate = 0.1,
                        depth = 3, iterations = 100, subsample = 1)
    )
  }

  checkmate::assert_list(gradient_params, names = "named")
  checkmate::assert_flag(na.rm)
  checkmate::assert_integerish(seed, lower = 1, null.ok = TRUE)
  checkmate::assert_integerish(bootstrap_cores, lower = 0)
  checkmate::assert_flag(mse)
  checkmate::assert_integerish(B, lower = 0)

  if (!is.null(seed)) set.seed(seed)

  if (!all(table(pop_data[, dom_name]) > 0))
    stop("All population domains must have at least one unit.")

  smp_check       <- smp_data[, setdiff(names(smp_data), dom_name), drop = FALSE]
  has_categorical <- any(sapply(smp_check, is.factor) | sapply(smp_check, is.character))
  if (has_categorical) {
    if (gbm_engine %in% c("xgboost", "lightgbm"))
      warning("One-hot encoding applied for categorical features (xgboost/lightgbm).")
    else if (gbm_engine == "catboost")
      message("CatBoost: native categorical handling enabled.")
  }

  list(smp_data = smp_data, pop_data = pop_data, gradient_params = gradient_params)
}


# Internal core MEGB fitter (renamed from MEGB::megb to avoid namespace collision)
#' @importFrom dplyr group_by summarise
#' @importFrom stats predict

megb_em <- function(Y, X, dom_name, smp_data, pop_data,
                    gradient_params  = NULL,
                    na.rm            = TRUE,
                    seed,
                    mse              = FALSE,
                    B                = 100,
                    bootstrap_cores  = 0,
                    gbm_engine       = "xgboost",
                    ...) {

  call       <- match.call()
  ts_gradient <- Sys.time()

  checked_inputs <- input_checks_megb(
    Y = Y, X = X, dom_name = dom_name, smp_data = smp_data,
    pop_data = pop_data, gradient_params = gradient_params,
    na.rm = na.rm, seed = seed, bootstrap_cores = bootstrap_cores,
    mse = mse, B = B, gbm_engine = gbm_engine
  )
  smp_data        <- checked_inputs$smp_data
  pop_data        <- checked_inputs$pop_data
  gradient_params <- checked_inputs$gradient_params
  rm(checked_inputs)

  if (na.rm) {
    comp_smp <- complete.cases(smp_data)
    smp_data <- smp_data[comp_smp, ]
    Y        <- Y[comp_smp]
    X        <- X[comp_smp, , drop = FALSE]
  }

  formula_random_effects <- paste0("(1|", dom_name, ")")
  cov_names              <- names(X)

  if (gbm_engine %in% c("xgboost", "lightgbm")) {
    # Build a model.frame ourselves with na.action = na.pass so NAs propagate
    # through one-hot encoding for xgboost/lightgbm to handle natively.
    # model.matrix.default does NOT forward na.action through ... to its
    # internal model.frame call — so passing na.action directly to
    # model.matrix() is silently ignored and rows still get dropped via the
    # global getOption("na.action"). Doing the model.frame step explicitly
    # produces an object with a "terms" attribute, which model.matrix then
    # uses as-is without re-calling model.frame.
    mm_formula  <- terms(~ . - 1, data = smp_data[, cov_names])
    mf_smp      <- stats::model.frame(mm_formula, data = smp_data,
                                       na.action = stats::na.pass)
    mf_pop      <- stats::model.frame(mm_formula, data = pop_data,
                                       na.action = stats::na.pass)
    X_smp_mm    <- stats::model.matrix(mm_formula, data = mf_smp)
    X_pop_mm    <- stats::model.matrix(mm_formula, data = mf_pop)
    X           <- X_smp_mm
    smp_data    <- cbind(smp_data[dom_name], as.data.frame(X_smp_mm))
    pop_data    <- cbind(pop_data[dom_name], as.data.frame(X_pop_mm))
    cov_names   <- colnames(X_smp_mm)
  }

  model <- em_gb_lmm(
    Y                      = Y,
    X                      = X,
    formula_random_effects = formula_random_effects,
    gradient_params        = gradient_params,
    data                   = smp_data,
    initial_random_effects = 0,
    max_iterations         = 25,
    error_tolerance        = 1e-04,
    dom_name               = dom_name,
    cov_names              = cov_names,
    gbm_engine             = gbm_engine,
    ...
  )

  unit_level_predictions <- gbm_predict(
    model      = model,
    smp_data   = smp_data,
    pop_data   = pop_data,
    Y          = Y,
    dom_name   = dom_name,
    gbm_engine = gbm_engine,
    cov_names  = cov_names
  )

  unit_pred_smp <- unit_level_predictions$unit_pred_smp
  gb_smp        <- unit_level_predictions$gb_smp
  unit_preds    <- unit_level_predictions$unit_pred_pop

  mean_preds <- unit_preds |>
    dplyr::group_by(dom_name) |>
    dplyr::summarise(Mean = mean(unit_preds)) |>
    as.data.frame()

  data_sum <- data_info(dom_name = dom_name, pop = pop_data, smp = smp_data)

  # ── n_d diagnostic: helps detect ML bias in bootstrap σ_u ──────────────────
  # Drop zero-count factor levels (out-of-sample domains) so statistics only
  # describe in-sample wards.
  n_d_tab    <- table(smp_data[[dom_name]])
  n_d_insamp <- as.numeric(n_d_tab[n_d_tab > 0])
  pct_sing   <- round(100 * mean(n_d_insamp == 1), 1)
  message(sprintf(
    "In-sample domains — n: %d  median n_d: %.1f  min: %d  max: %d  n_d=1: %d (%.1f%%)",
    length(n_d_insamp), median(n_d_insamp),
    min(n_d_insamp), max(n_d_insamp),
    sum(n_d_insamp == 1), pct_sing
  ))
  message(sprintf(
    "LMM variance components — σ_u: %.4f  σ_e: %.4f",
    model$ran_eff_sd, model$error_sd
  ))

  if (mse) {
    message("Bootstrap with ", B, " iterations has started")
    mse_estimated <- mse_megb(
      Y                      = Y,
      X                      = X,
      dom_name               = dom_name,
      smp_data               = smp_data,
      model                  = model,
      error_sd               = model$error_sd,
      pop_data               = pop_data,
      B                      = B,
      initial_random_effects = 0,
      ErrorTolerance         = 0.0001,
      MaxIterations          = 10,
      cov_names              = cov_names,
      gradient_params        = gradient_params,
      formula_random_effects = formula_random_effects,
      bootstrap_cores        = bootstrap_cores,
      seed                   = seed,
      gbm_engine             = gbm_engine,
      unit_pred_smp          = unit_pred_smp,
      gb_smp                 = gb_smp,
      unit_preds             = unit_preds,
      ...
    )
  } else {
    mse_estimated <- NULL
  }

  res <- list(
    call          = call,
    Indicators    = mean_preds,
    data_sum      = data_sum,
    megb_model    = model,
    inp_smp_data  = list(smp_data = smp_data, target_var = Y),
    unit_preds_all = unit_preds,
    unit_pred_smp  = unit_pred_smp,
    gb_smp         = gb_smp,
    MSE_Estimates  = mse_estimated,
    time_gradient  = Sys.time() - ts_gradient,
    gradient_params = gradient_params,
    gbm_engine      = gbm_engine,
    X_proc                 = X,
    pop_data_proc          = pop_data,
    cov_names_proc         = cov_names,
    formula_random_effects = formula_random_effects
  )
  class(res) <- "MEGB"
  res
}
