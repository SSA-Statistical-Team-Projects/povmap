
#' MEGB
#'
#' Gradient Boosting parameter,
#' @return list of class MEGB with results

# Inputs: Gradient Boosting parameter
params <- list(
  eta = 0.1,              # learning rate
  nrounds = 3000,          # max number of trees
  max_depth = 5,           # max depth of trees
  min_child_weight = 5,    # minimal weight sum for child nodes
  subsample = 0.8,         # subsample ratio per tree
  colsample_bytree = 0.8,  # subsample ratio of columns per tree
  lambda = 0,              # L2 regularisation
  alpha = 0.1              # L1 regularisation
)

#' This function defines input checks for Mixed Effect Gradient Boosting,
#' which combines Gradient Boosting with Linear Mixed Models (LMM).
#' @param Y numeric, target variable
#' @param X data.frame, Covariates
#' @param dom_name factor, name of the domain variable
#' @param smp_data data.frame, sample data
#' @param pop_data data.frame, population data
#' @param gradient_params list, optional list of hyperparameters passed to the selected GB engine and required to match the specification of gbm_engine.
#' @param na.rm logical, should NAs be removed?
#' @param seed integer, seed to reproduce results
#' @param bootstrap_cores integer, number of cores to use for bootstrap
#' @param mse logical, should MSE be estimated?
#' @param B integer, number of bootstrap samples for MSE
#' @param gbm_engine character, choice of Gradient Boosting implementation ("xgboost", "lightgbm", "catboost")
#' @return logical to indicate if tests were passed or not
#' @import checkmate


input_checks_megb <- function(Y, X, dom_name, smp_data, pop_data, 
                              gradient_params, na.rm, seed, bootstrap_cores, 
                              mse, B, gbm_engine) {

  # Check Y
  assert_numeric(Y, any.missing = FALSE, min.len = 1, finite = TRUE, 
                 null.ok = FALSE, len = nrow(X))

  # Check X
  assert_data_frame(X, col.names =  "named", min.rows = 1, min.cols = 1, 
                    any.missing = FALSE)
  # Ensure X has column names
  assert_character(names(X), any.missing = FALSE, min.len = 1)

  # Check smp_data and pop data
  assert_data_frame(pop_data, min.rows = 1, any.missing = TRUE)
  assert_data_frame(smp_data, min.rows = 1, any.missing = TRUE)
  assert_names(intersect(colnames(pop_data), colnames(smp_data)), must.include = c(dom_name))


    # check character vars in smp_data
    for (col_name in names(smp_data)) {
    # check if variable has class "character"
    if (is.character(smp_data[[col_name]])) {
      # transform to factor
      smp_data[[col_name]] <- as.factor(smp_data[[col_name]])
      warning("The variable ", col_name, " of smp_data was transformed from character to factor.")
    }
    }

    # check domain-variable is factor
  assert_factor(smp_data[[dom_name]], any.missing = FALSE)


    # check character vars in pop_data
    for (col_name in names(pop_data)) {
    # check if variable has class "character"
    if (is.character(pop_data[[col_name]])) {
      # transform to factor
      pop_data[[col_name]] <- as.factor(pop_data[[col_name]])
      warning("The variable ", col_name, " of pop_data was transformed from character to factor.")
    }
    }

    # check domain-variable is factor
  assert_factor(pop_data[[dom_name]], any.missing = FALSE)

  # Ensure all variables in smp_data exist in pop_data (except for Y)
  Y_name <- tail(strsplit(deparse(substitute(Y)), "\\$")[[1]], 1)
  smp_vars <- colnames(smp_data)[colnames(smp_data) != Y_name]
  missing_vars <- setdiff(smp_vars, colnames(pop_data))
  if (length(missing_vars) > 0) {
    stop(paste("The following variables from smp_data are missing in pop_data:",
               paste(missing_vars, collapse = ", ")))
  }

  # Ensure there are no domains in smp_data that are not in pop_data
  if (!all(unique(smp_data[[dom_name]]) %in% unique(pop_data[[dom_name]]))) {
    stop("There are domains in smp_data that are not present in pop_data.")
  }

  # check gbm engine is indicated
  assert_choice(gbm_engine, choices = c("xgboost", "lightgbm", "catboost"))
  
  # check if catboost package is required and installed
  if (gbm_engine == "catboost" && !requireNamespace("catboost", quietly = TRUE)) {
    stop("The 'catboost' package is required for gbm_engine = 'catboost'. Please
         install it.")
  }
  
  # set default parameter depending on gbm_engine
  if (is.null(gradient_params)) {
    
    message("gradient_params is NULL — default parameters are used for the specified gbm_engine.")
    
    gradient_params <- switch(gbm_engine,
                    "xgboost"  = list(
                      eta = 0.1,
                      max_depth = 3,
                      nrounds = 100,
                      subsample = 1
                    ),
                    "lightgbm" = list(
                      objective = "regression",
                      metric = "rmse", 
                      learning_rate = 0.1,
                      num_leaves = 63,
                      n_estimators = 100,
                      bagging_fraction = 1,
                      nrounds = 100
                    ),
                    "catboost" = list(
                      loss_function = "RMSE",
                      learning_rate = 0.1,
                      depth = 3,
                      iterations = 100,
                      subsample = 1
                    )
    )
  }
  
  # Check gradient_params
  assert_list(gradient_params, names = "named")

  # Check na.rm
  assert_flag(na.rm)

  # Check seed
  assert_integerish(seed, lower = 1, null.ok = TRUE)

  # Check cores
  assert_integerish(bootstrap_cores, lower = 0)

  # Check mse
  assert_flag(mse)

  # Check B
  assert_integerish(B, lower = 0)

  # If seed is provided, set it for reproducibility
  if (!is.null(seed)) {
    set.seed(seed)
  }
  if(!all(table(pop_data[, dom_name]) > 0)){
    stop("Domains that are not in population data can not be predicted.")
  }

  # Check for categorial variables 
  smp_check <- smp_data[ , setdiff(names(smp_data), dom_name), drop = FALSE]
  
  has_categorical <- any(sapply(smp_check, is.factor) | sapply(smp_check, is.character))
  
  if (has_categorical) {
    if (gbm_engine %in% c("xgboost", "lightgbm")) {
      warning("Using XGBoost/LightGBM: applying one-hot encoding for categorical features")
    } else if (gbm_engine == "catboost") {
      message("Using CatBoost: native categorical handling enabled (no one-hot).")
    }
  }

  
  # return a list with changed smp_data and pop_data and gradient_params  
  list(
    smp_data = smp_data,
    pop_data = pop_data,
    gradient_params = gradient_params
  )
  

}

#' This function defines Mixed Effect Gradient Boosting,
#' which combines Gradient Boosting with Linear Mixed Models (LMM).
#' @param Y numeric, target variable
#' @param X data.frame, Covariates
#' @param dom_name factor, name of the domain variable
#' @param smp_data data.frame, sample data
#' @param pop_data data.frame, population data
#' @param gradient_params list, optional list to define Gradient Boosting parameters. If used, has to use the parameter names of the used GB engine.
#' @param na.rm logical, should NAs be removed?
#' @param seed integer, seed to reproduce results
#' @param mse logical, should MSE be estimated?
#' @param B numeric, number of bootstrap samples for MSE
#' @param bootstrap_cores integer, number of cores to use
#' @param gbm_engine character, choice of Gradient Boosting implementation ("xgboost", "lightgbm", "catboost")
#' @param ... additional arguments
#' @return list of class MEGB with results:
#' \itemize{
#'        \item Indicators, containing the average domain predictions
#'        \item gbmodel, containing the model details
#'        \item inp_smp_data, containing the input data
#'        \item unit_preds_all, containing the unit-level predictions (population)
#'        \item unit_preds_all, containing the unit-level predictions (sample)
#'        \item MSE_Estimates, containing te estimates MSE values
#'        \item time_gradient, containing the runtime
#'        \item gradient_params, containing the used Gradient Boosting parameters
#' }
#' @examples
#' \donttest{
#' # Loading data - population and sample data
#' data("population_data")
#' data("sample_data")
#'
#' # Example 1: with MSE = TRUE
#'
#' params <- list(eta = 0.1,      # learning rate
#'                nrounds = 100   # number of iterations
#'                )
#'
#' result <- megb(Y = sample_data$y,
#'                X = sample_data[, 2:6],
#'                dom_name = "domain",
#'                smp_data = sample_data,
#'                pop_data = population_data,
#'                gradient_params = params,
#'                na.rm = TRUE,
#'                seed = 1,
#'                mse = TRUE,
#'                B = 2,
#'                bootstrap_cores = 0,
#'                gbm_engine = "xgboost"
#'                )
#'}
#' @import lme4
#' @importFrom xgboost xgboost xgb.cv xgb.importance
#' @importFrom dplyr group_by summarise `%>%`
#' @export

# main function: Mixed Effect Gradient Boosting
megb <- function(Y,                  # target variable
                 X,                  # Covariates (dataframe)
                 dom_name,           # name of domain variable (e.g. groups)
                 smp_data,           # sample data
                 pop_data,           # population data (for predictions)
                 gradient_params = NULL,  # optional: Gradient Boosting parameters
                 na.rm = TRUE,       # should NAs be removed?
                 seed,               # seed to reproduce results
                 mse = FALSE,        # should MSE be estimated?
                 B = 100,            # number of bootstrap samples for MSE
                 bootstrap_cores = 0,# number of cores
                 gbm_engine = "xgboost", # choice of Gradient Boosting implementation
                 ...){               # additional arguments

  call <- match.call()
  # timestamp to measure runtime
  ts_gradient <- Sys.time()

  checked_inputs <- input_checks_megb(Y = Y,
                    X = X,
                    dom_name = dom_name,
                    smp_data = smp_data,
                    pop_data = pop_data,
                    gradient_params = gradient_params,
                    na.rm = na.rm,
                    seed = seed,
                    bootstrap_cores = bootstrap_cores,
                    mse = mse,
                    B = B,
                    gbm_engine = gbm_engine
                    )

# Save the changed inputs
  smp_data <- checked_inputs$smp_data
  pop_data <- checked_inputs$pop_data
  gradient_params <- checked_inputs$gradient_params
  rm(checked_inputs)

  # optional: remove NAs
  if (na.rm) {
    comp_smp <- complete.cases(smp_data)
    smp_data <- smp_data[comp_smp, ]
    Y <- Y[comp_smp]
    X <- X[comp_smp, , drop = FALSE]
  }

  # create formula for Random Effects
  formula_random_effects <- paste0("(1|", dom_name, ")")
  cov_names <- names(X)

  # --------------------------------------------------------------
  # One-hot encoding (only for XGBoost / LightGBM)
  # --------------------------------------------------------------
  if (gbm_engine %in% c("xgboost", "lightgbm")) {

    mm_formula <- terms( ~ . - 1, data = smp_data[, cov_names])
    # Design-matrices with identical structure
    X_smp_mm <- model.matrix(mm_formula, data = smp_data)
    X_pop_mm <- model.matrix(mm_formula, data = pop_data)
    
    # Update for X, smp_data, pop_data
    X        <- X_smp_mm
    smp_data <- cbind(smp_data[dom_name], as.data.frame(X_smp_mm))
    pop_data <- cbind(pop_data[dom_name], as.data.frame(X_pop_mm))
    
    # update covariable names
    cov_names <- colnames(X_smp_mm)
  }
  
  
  # calculate point-estimates ---------------------------------------------
  # EM-algorithm: Combination of Gradient Boosting and LMM
  model <- em_gb_lmm(
    Y = Y,
    X = X,
    formula_random_effects = formula_random_effects,
    gradient_params = gradient_params,
    data = smp_data,
    initial_random_effects = 0, # initial values for domain-effects
    max_iterations = 10,        # max EM-iterations
    error_tolerance = 1e-04,    
    dom_name = dom_name,
    cov_names = cov_names,
    gbm_engine = gbm_engine,
    ...
  )
  unit_level_predictions <- gbm_predict(
    model = model,
    smp_data = smp_data,
    pop_data = pop_data,
    Y = Y,
    dom_name = dom_name,
    gbm_engine = gbm_engine,
    cov_names = cov_names
  )

  # predictions for sample data
  unit_pred_smp <- unit_level_predictions$unit_pred_smp
  # residuals for sample data
  res <- unit_level_predictions$res
  # predictions for pop data
  unit_preds <- unit_level_predictions$unit_pred_pop

  # average prediction for each domain
  mean_preds <- unit_preds %>%
    group_by(dom_name) %>%
    dplyr::summarise(Mean = mean(unit_preds))
  
  data_sum <- data_info(dom_name = dom_name, pop = pop_data, smp = smp_data)

  # MSE-Predictions ---------------------------------------------------------
  if (mse) {
    message(paste("Bootstrap with", B, "iterations has started"))
    mse_estimated <- mse_megb(
      Y = Y,
      X = X,
      dom_name = dom_name,
      smp_data = smp_data,
      model = model,
      error_sd = model$error_sd,
      pop_data = pop_data,
      B = B,
      initial_random_effects = 0,    # initial value for Random Effects
      ErrorTolerance = 0.0001,     # error tolerance for EM-algorithm
      MaxIterations = 10,
      cov_names = cov_names,
      gradient_params = gradient_params,
      formula_random_effects = formula_random_effects,
      bootstrap_cores = bootstrap_cores,
      seed = seed, 
      gbm_engine = gbm_engine,
      unit_pred_smp = unit_pred_smp,
      unit_preds = unit_preds, 
      ...
    )
  } else {
    mse_estimated <- NULL
  }

  # list of results
  res <- list(
    call = call,
    Indicators = mean_preds,          # avg predictions per domain
    data_sum = data_sum,
    megb_model = model, # model details
    inp_smp_data = list(smp_data = smp_data, target_var = Y), # input_data
    unit_preds_all = unit_preds,      # unit-level predictions (population)
    unit_pred_smp = unit_pred_smp,    # unit-level predictions (sample)
    MSE_Estimates = mse_estimated,    # estimated MSE values
    time_gradient = Sys.time() - ts_gradient, # runtime
    gradient_params = gradient_params, # used Gradient Boosting parameters
    gbm_engine = gbm_engine
  )

  # Class output as "MEGB"
  class(res) <- "MEGB"
  res
}


