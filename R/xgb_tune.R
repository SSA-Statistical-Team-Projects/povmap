#' Tuning extreme gradient boosting models for domain-level averages
#'
#' The funtion \code{xgb_tune} fine-tunes the hyperparameters for extreme gradient boosting models, following \cite{Merfeld and Newhouse (2023)}.
#' It offers the flexibility to allocate \code{domains} to folds and utilizes estimated
#' means at the domain-level for cross-validation. Users can specify the number of
#' folds, including an option for leave-one-out cross-validation.
#'
#' @param fixed a two-sided linear formula object describing the
#' fixed-effects part of the model with the dependent variable on the left
#' of a ~ operator and the explanatory variables on the right, separated
#' by + operators. All variables (except for \code{domains} and \code{cluster})
#' must be numeric.
#' @param smp_data a data frame that needs to comprise all variables including
#' \code{domains} and \code{cluster}.
#' @param smp_weights a character string containing the name of the variable that
#' indicates weights in \code{smp_data}. The variable has to be numeric.
#' Defaults to \code{NULL}.
#' @param domains a character string containing the name of a variable
#' that indicates domains in \code{smp_data}. The variable can be
#' numeric or a factor.
#' @param cluster a character string containing the name of a variable
#' that indicates clusters in \code{smp_data}. The variable can be
#' numeric or a factor. Defaults to \code{"domains"}.
#' @param transformation a character string. Two different transformation
#' types for the dependent variable can be chosen (i) no transformation ("no");
#' (ii) log transformation ("log"); (iii) Arcsin transformation ("arcsin").
#' Defaults to \code{"no"}.
#' @param folds number of folds. Defaults to 10.
#' @param nround combination of maximum number of boosting iterations. Defaults to 150 and 250.
#' @param max_depth combination of maximum depth of a tree. Defaults to 4 and 6.
#' @param colsample_bytree combination of subsample ratios of columns when constructing each tree. Defaults to 0.6 and 1.
#' @param colsample_bylevel combination of subsample ratios of columns for each level. Defaults to 0.6 and 1.
#' @param colsample_bynode combination of subsample ratios of columns for each node (split). Defaults to 0.6 and 1.
#' @param subsample combination of subsample ratios of the training instances. Defaults to 0.6 and 1.
#' @param min_child_weight minimum sum of instance weight required in a child node.
#' If the tree partitioning step produces a leaf node with a sum of instance weight
#' less than \code{min_child_weight}, then the building process will cease further
#' partitioning. A larger value of \code{min_child_weight} leads to a more conservative
#' algorithm. Defaults to 1.
#' @param eta step size shrinkage. After each boosting step, one can obtain the
#' weights of new features directly, and the parameter \code{eta} is used to shrink
#' these feature weights, thereby making the boosting process more conservative.
#' Range of [0, 1]. Defaults to 0.3.
#' @param gamma minimum loss reduction needed to create an additional partition
#' on a leaf node of the tree. A larger value of \code{gamma} corresponds to a more
#' conservative algorithm. Defaults to 0.
#' @param max_delta_step maximum allowed step size for adjusting the output of each
#' leaf. If the value is set to 0, it indicates that there is no constraint Defaults to 0.
#' @param lambda L2 regularization term on weights. Increasing this value will result in a more conservative model.
#' Defaults to 1.
#' @param alpha L1 regularization term on weights. Increasing this value will result in a more conservative model.
#' @param search search strategy. \code{"grid"} (the default) evaluates the
#'   full Cartesian product of the parameter vectors above and reproduces the
#'   behaviour of earlier versions exactly. \code{"random"} instead draws
#'   \code{n_iter} configurations from continuous ranges given in \code{bounds}.
#' @section Choosing a search strategy:
#' Random search was measured against the grid and against Bayesian
#' optimisation on a sub-area estimation problem with 10,979 sub-areas, 636
#' municipalities and 23 covariates, scoring domain-level weighted MSE under
#' grouped CV.
#'
#' Random search at 63 s beat the 32-point grid at 61 s by 23 percent on
#' out-of-sample R2 (0.139 against 0.113), because the grid tests two values of
#' five parameters and so establishes direction rather than locating an
#' optimum. Adding \code{nround} to the search improved R2 by a further 10
#' percent (0.153) but cost four times the wall clock, so the cheap
#' configuration reaches within 9 percent of the expensive one. If you are
#' wall-clock constrained, the first minute buys most of the available gain.
#'
#' Bayesian optimisation (ParBayesianOptimization) was tested and LOST to
#' random search at equal wall clock, which is why no such option is offered
#' here. A CV evaluation on this problem costs about 6 s while the surrogate
#' costs 25 to 29 s per iteration, so the machinery is four to five times more
#' expensive than the work it economises on; at equal wall clock random search
#' affords roughly six times as many evaluations. The surrogate also failed an
#' informativeness check, twenty guided iterations improving on ten random
#' initial points by 0.3 percent and by nothing on a second seed. The
#' conclusion is specific to cheap evaluations and would likely invert if a
#' single evaluation took minutes.
#' @param n_iter number of configurations drawn when \code{search = "random"}.
#'   Defaults to 32. Note that 32 random draws is NOT equivalent in coverage to a
#'   32 point grid: the grid tests two values of five parameters and so
#'   establishes direction rather than locating an optimum, whereas the random
#'   draws explore the interior of the ranges but cover any single axis more
#'   thinly. Raising \code{n_iter} trades compute for resolution directly.
#' @param bounds named list giving the search range per parameter, used only when
#'   \code{search = "random"}. Each element is either a bare \code{c(lo, hi)},
#'   which inherits that parameter's default sampling scale, or
#'   \code{list(range = c(lo, hi), scale = "log")} to override it. Only
#'   parameters named here are varied; every other parameter is held at the first
#'   value of its own argument, so exploring depth does not silently move
#'   \code{gamma} and \code{alpha} as well. Defaults to \code{eta},
#'   \code{max_depth}, \code{min_child_weight}, \code{subsample},
#'   \code{colsample_bytree} and \code{lambda}. \code{eta},
#'   \code{min_child_weight} and \code{lambda} are drawn log-uniformly because
#'   they span orders of magnitude and a uniform draw would spend most of its
#'   budget at the insensitive end; \code{alpha} is drawn with a point mass at
#'   zero, since a log scale is impossible at zero.
#' @param early_stopping_rounds if supplied, each configuration is fitted with
#'   \code{nrounds_max} rounds and stopped early on an inner validation group
#'   held out of the training folds and grouped by \code{cluster}. This takes
#'   \code{nround} out of the search entirely, freeing a dimension so a fixed
#'   budget goes further. The inner group is separate from the outer fold, which
#'   remains purely for scoring. Defaults to \code{NULL}, i.e. fixed
#'   \code{nround} as before.
#' @param nrounds_max ceiling on boosting rounds when early stopping is in use.
#'   Defaults to 2000. If any fold stops at this ceiling a warning is raised,
#'   since the budget rather than the data then chose the stopping point.
#' @param seed integer seed making the whole routine reproducible: the assignment
#'   of clusters to folds, the random draw of configurations, and the underlying
#'   \code{xgboost} fits, which are otherwise stochastic whenever
#'   \code{subsample} or the \code{colsample_*} parameters are below 1. The fit
#'   seed is passed straight to xgboost, so configurations share common random
#'   numbers. The fold assignment and the xgboost fits were both previously
#'   unseeded, so results varied run to run even under grid search;
#'   \code{seed = NULL} preserves that behaviour exactly.
#' Defaults to 0.
#' @param cpus. Number of cores to parallelize across. Defaults to 1 (no parallelization)
#' @param verbose display progress. Defaults to FALSE.
#' @param ... additional parameters to be passed to \code{xgb.train}.
#'
#' @return An object of class \code{xgb}, \code{emdi}, containing the optimal
#' hyperparameters for an extreme gradient boosting model and the out-of-sample
#' R-squared at the domain level (\code{r2_oos}).
#' @references
#' Merfeld, J. D., & Newhouse, D. (2023). Improving Estimates of Mean Welfare and Uncertainty
#' in Developing Countries (No. 10348). The World Bank. \cr \cr
#' @export
#' @importFrom xgboost xgboost
#' @importFrom dplyr left_join
#' @importFrom foreach foreach %dopar%
#' @importFrom doParallel registerDoParallel stopImplicitCluster
#'
#' @examples
#' \donttest{
#' # Loading data - population and sample data
#' data("eusilcA_pop")
#' data("eusilcA_smp")
#'
#' xgb_tune_model <- xgb_tune(fixed = eqIncome ~ eqsize + cash + self_empl +
#'                            unempl_ben + age_ben + surv_ben + sick_ben +
#'                            dis_ben + rent + fam_allow + house_allow +
#'                            cap_inv + tax_adj + district,
#'                            smp_data = eusilcA_smp,
#'                            domains = "district")
#'}

xgb_tune <- function(fixed,
                     smp_data,
                     smp_weights = NULL,
                     domains,
                     cluster = "domains",
                     transformation = "no",
                     folds = 10,
                     nround = c(150, 300),
                     max_depth = c(3, 4),
                     colsample_bytree = c(0.6, 0.8),
                     colsample_bylevel = c(1),
                     colsample_bynode = c(1),
                     subsample = c(0.8),
                     min_child_weight = c(5, 15),
                     eta = c(0.3),
                     gamma = c(0),
                     max_delta_step = c(0),
                     lambda = c(0.5, 1.5),
                     alpha = c(0),
                     search = c("grid", "random"),
                     n_iter = 32,
                     bounds = NULL,
                     early_stopping_rounds = NULL,
                     nrounds_max = 2000,
                     seed = NULL,
                     cpus = 1, 
                     verbose = TRUE,
                     ...){

  # Pinned-version guard (see BUILD_PIN_xgb.txt) -- fail before any model fit
  .assert_xgb_version()

  # Reject anything that fell into `...` unused. `...` is forwarded to
  # xgb.train(), so a caller written against a newer signature -- or a
  # misspelled formal -- would otherwise be discarded in silence and the search
  # would quietly run on defaults. See .check_xgb_dots() for the full reasoning.
  .check_xgb_dots(list(...), fn_formals = names(formals(sys.function())))

  # Data preparation
  #_____________________________________________________________________________
  outcome <- all.vars(fixed[[2]])
  covariates <- all.vars(fixed[[3]])

  #split <- strsplit(as.character(fixed), "~", fixed = TRUE)
  #outcome <- trimws(split[[1]][1])
  #covariates <- trimws(strsplit(trimws(split[[1]][2]), "\\+")[[1]])
  X_smp <- smp_data[,c(covariates,domains)]
  Y_smp <- data.frame(smp_data[,outcome])

  if(is.null(smp_weights)==FALSE){
    smp_weights <- smp_data[,smp_weights]
  } else {
    smp_weights <- rep(1, length = nrow(Y_smp))
  }
  if (cluster=="domains"){
    cluster <- paste0(domains)
  }
  colnames(Y_smp) <- "labels"

  # Check
  #_____________________________________________________________________________
  xgb_check2(
    transformation = transformation,
    Y_smp = Y_smp,
    X_smp = X_smp,
    smp_weights = smp_weights,
    domains = domains,
    cluster = cluster)

  # Transformation
  #_____________________________________________________________________________
  if (transformation=="arcsin"){
    Y_smp <- asin(sqrt(Y_smp))
  }
  if (transformation=="log"){
    Y_smp <- log(Y_smp)
  }

  # Search mode, seed
  #_____________________________________________________________________________
  search <- match.arg(search)
  ## Seeding covers BOTH the fold assignment below and the random draw. The fold
  ## assignment was previously unseeded, so even grid search was not reproducible
  ## run to run; seed = NULL preserves that old behaviour exactly.
  if (!is.null(seed)) set.seed(seed)

  # Folds
  #_____________________________________________________________________________
  cluster_col <- data.frame(X_smp[[paste0(cluster)]])
  colnames(cluster_col) <- paste0(cluster)
  cluster_unique <- data.frame(unique(cluster_col[,1]))
  colnames(cluster_unique) <- paste0(cluster)
  cluster_unique$fold <- sample(x = 1:folds, size = nrow(cluster_unique), replace = TRUE)

  cluster_col <- cluster_col %>%
    dplyr::left_join(cluster_unique, by = paste0(cluster))
  if (cluster=="domains"){
    X_final <- X_smp[,-c(which(colnames(X_smp)==paste0(cluster)))]
  } else{
    X_final <- X_smp[,-c(which(colnames(X_smp)==paste0(cluster)), which(colnames(X_smp)==paste0(domains)))]
  }

  # Grid
  #_____________________________________________________________________________
  if (search == "grid") {
  tunegrid <- expand.grid(
    nround             = nround,
    max_depth          = max_depth,
    colsample_bytree   = colsample_bytree,
    colsample_bylevel  = colsample_bylevel,
    colsample_bynode   = colsample_bynode,
    subsample          = subsample,
    min_child_weight   = min_child_weight,
    eta                = eta,
    gamma              = gamma,
    max_delta_step     = max_delta_step,
    lambda             = lambda,
    alpha              = alpha
  )
  } else {
    ## Random search: sample the named parameters from continuous ranges and hold
    ## every other parameter at the FIRST value of its argument.
    if (is.null(bounds)) bounds <- .XGB_DEFAULT_BOUNDS
    if (!is.list(bounds) || is.null(names(bounds)) || any(!nzchar(names(bounds))))
      stop("bounds must be a NAMED list, one entry per parameter to vary.", call. = FALSE)
    known <- c("nround", "max_depth", "colsample_bytree", "colsample_bylevel",
               "colsample_bynode", "subsample", "min_child_weight", "eta",
               "gamma", "max_delta_step", "lambda", "alpha")
    bad <- setdiff(names(bounds), known)
    if (length(bad))
      stop("bounds names not recognised: ", paste(bad, collapse = ", "), call. = FALSE)
    sampled <- .xgb_sample_bounds(bounds, n_iter)
    fixed_vals <- list(
      nround = nround[1], max_depth = max_depth[1],
      colsample_bytree = colsample_bytree[1], colsample_bylevel = colsample_bylevel[1],
      colsample_bynode = colsample_bynode[1], subsample = subsample[1],
      min_child_weight = min_child_weight[1], eta = eta[1], gamma = gamma[1],
      max_delta_step = max_delta_step[1], lambda = lambda[1], alpha = alpha[1])
    tunegrid <- as.data.frame(lapply(known, function(p)
      if (p %in% names(sampled)) sampled[[p]] else rep(fixed_vals[[p]], n_iter)))
    names(tunegrid) <- known
  }

  ## Parallel backend, shared by every search mode.
  cl <- parallel::makeCluster(cpus)
  doParallel::registerDoParallel(cl)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  dots <- list(...)

  score_args <- list(X_final = X_final, X_smp = X_smp, Y_smp = Y_smp,
                     smp_weights = smp_weights, cluster_col = cluster_col,
                     cluster = cluster, domains = domains, folds = folds,
                     early_stopping_rounds = early_stopping_rounds,
                     nrounds_max = nrounds_max, seed = seed, dots = dots)

  sc <- do.call(.xgb_score_configs, c(list(tunegrid = tunegrid), score_args))
  OPT <- sc$OPT; ITER <- sc$ITER; domain_labels_list <- sc$domain_labels

  # Optimal values
  #_____________________________________________________________________________
  mean_mse <- apply(OPT, 1, FUN = mean)
  best_row <- which.min(mean_mse)
  mse_min  <- mean_mse[best_row]

  es_on <- !is.null(early_stopping_rounds)
  ## With early stopping the useful nrounds is the one each fold learned, not
  ## the grid value, so report the mean stopping point for the selected config.
  best_iters <- ITER[best_row, ]
  nround_opt <- if (es_on) round(mean(best_iters)) else tunegrid$nround[best_row]
  max_depth_opt <- tunegrid$max_depth[best_row]
  colsample_bytree_opt <- tunegrid$colsample_bytree[best_row]
  colsample_bylevel_opt <- tunegrid$colsample_bylevel[best_row]
  colsample_bynode_opt <- tunegrid$colsample_bynode[best_row]
  subsample_opt <- tunegrid$subsample[best_row]
  min_child_weight_opt <- tunegrid$min_child_weight[best_row]
  eta_opt <- tunegrid$eta[best_row]
  max_delta_step_opt <- tunegrid$max_delta_step[best_row]
  gamma_opt <- tunegrid$gamma[best_row]
  lambda_opt <- tunegrid$lambda[best_row]
  alpha_opt <- tunegrid$alpha[best_row]

  # Out-of-sample R2 at the domain level
  # Variance of domain-level mean outcomes pooled across all folds
  all_domain_labels <- unlist(domain_labels_list)
  var_labels <- var(all_domain_labels)
  r2_oos <- 1 - mse_min / var_labels

  final_output <- list(nround_opt, max_depth_opt, colsample_bytree_opt, colsample_bylevel_opt,
                       colsample_bynode_opt, subsample_opt, min_child_weight_opt, eta_opt,
                       gamma_opt, max_delta_step_opt, lambda_opt, alpha_opt, mse_min, r2_oos)
  names(final_output) <- c("nround", "max_depth", "colsample_bytree", "colsample_bylevel",
                           'colsample_bynode', "subsample", "min_child_weight", "eta",
                           "gamma", "max_delta_step", 'lambda', "alpha", "mse_oos", "r2_oos")

  ## Extra diagnostics, appended so positional access to the original 14
  ## elements is unaffected.
  final_output$search    <- search
  final_output$n_configs <- nrow(tunegrid)
  if (es_on) {
    final_output$best_iters_by_fold <- best_iters
    ## Hitting the ceiling is the early-stopping analogue of an optimum sitting
    ## on a grid edge: the budget bound, not the data, chose the stopping point.
    final_output$hit_nrounds_max <- any(best_iters >= nrounds_max)
    if (isTRUE(final_output$hit_nrounds_max))
      warning("At least one fold stopped at nrounds_max (", nrounds_max,
              "); raise nrounds_max, the stopping point was bounded by the budget.",
              call. = FALSE)
  }
  if (search == "random") {
    final_output$bounds   <- bounds
    final_output$n_iter   <- n_iter
    final_output$position <- .xgb_bounds_position(bounds, final_output)
    if (any(final_output$position$at_boundary))
      warning("Selected value sits on a boundary for: ",
              paste(final_output$position$parameter[final_output$position$at_boundary],
                    collapse = ", "),
              ". The range was too narrow -- widen it and rerun.", call. = FALSE)
  }

  class(final_output) <- c("xgb","emdi")
  return(final_output)
}



#' Tuning extreme gradient boosting models for domain-level averages
#'
#' The funtion \code{xgb_tune} fine-tunes the hyperparameters for extreme gradient boosting models, following \cite{Merfeld and Newhouse (2023)}.
#' It offers the flexibility to allocate \code{domains} to folds and utilizes estimated
#' means at the domain-level for cross-validation. Users can specify the number of
#' folds, including an option for leave-one-out cross-validation.
#'
#' @param fixed a two-sided linear formula object describing the
#' fixed-effects part of the model with the dependent variable on the left
#' of a ~ operator and the explanatory variables on the right, separated
#' by + operators. All variables (except for \code{domains} and \code{cluster})
#' must be numeric.
#' @param smp_data a data frame that needs to comprise all variables including
#' \code{domains} and \code{cluster}.
#' @param smp_weights a character string containing the name of the variable that
#' indicates weights in \code{smp_data}. The variable has to be numeric.
#' Defaults to \code{NULL}.
#' @param domains a character string containing the name of a variable
#' that indicates domains in \code{smp_data}. The variable can be
#' numeric or a factor.
#' @param cluster a character string containing the name of a variable
#' that indicates clusters in \code{smp_data}. The variable can be
#' numeric or a factor. Defaults to \code{"domains"}.
#' @param transformation a character string. Two different transformation
#' types for the dependent variable can be chosen (i) no transformation ("no");
#' (ii) log transformation ("log"); (iii) Arcsin transformation ("arcsin").
#' Defaults to \code{"no"}.
#'
#' @import dplyr rlang recipes


tidy_xgb_tune <- function(fixed,
                          smp_data,
                          smp_weights = NULL,
                          domains,
                          cluster = "domains",
                          transformation = "no",
                          folds = 10,
                          tune_params = list("trees" = c(1L, 2000L),
                                             "mtry" = NULL,
                                             "min_n" = c(2L, 40L),
                                             "tree_depth" = c(1L, 15L),
                                             "learn_rate" = c(-10, -1),
                                             "loss_reduction" = c(-10, 1.5),
                                             "sample_size" = c(0.1, 1),
                                             "stop_iter" = NULL,
                                             "alpha" = c(-10, 1),
                                             "lambda" = c(-10 ,1),
                                             "colsample_bytree" = c(0.5, 1),
                                             "colsample_bylevel" = c(0.5, 1),
                                             "colsample_bynode" = c(0.5, 1),
                                             "max_delta_step" = c(0, 10)),
                          tune_size = NULL,
                          verbose = TRUE,
                          parallel_over = NULL) {

  ## ---- data preparation ---- ##
  yvar_chr <- all.vars(fixed[[2]])
  xvar_list <- all.vars(fixed[[3]])


  # Create outcome and predictor objects
  xsmp <- smp_data[, xvar_list, drop = FALSE]
  ysmp <- smp_data[, yvar_chr, drop = FALSE]
  colnames(ysmp) <- "labels"

  # Handle weights
  if (!is.null(smp_weights)) {
    smp_data[["smp_weights"]] <- smp_data[[smp_weights]]
  } else {
    smp_data[["smp_weights"]] <- rep(1, nrow(ysmp))
  }

  smp_data[["domains"]] <- smp_data[[domains]]
  smp_data[["cluster"]] <- smp_data[[cluster]]

  ### drop the weights and domain variable
  smp_data[[domains]] <- NULL
  smp_data[[smp_weights]] <- NULL
  smp_data[[cluster]] <- NULL

  ## quickly update formula for later
  fixed <-
    paste(fixed,
          paste("smp_weights + domains"),
          sep = " + ") |>
    as.formula()


  sample_tbl <-
    smp_data |>
    as_tibble()

  ## ---- transformation ---- ##

  if (transformation == "arcsin") {

    sample_tbl <-
      sample_tbl |>
      mutate(!!yvar_chr := asin(sqrt(!!sym(yvar_chr))))

  } else if (transformation == "log") {

    sample_tbl <-
      sample_tbl |>
      mutate(!!yvar_chr := log(!!sym(yvar_chr)))

  }


  ## ---- preprocessing and data preparation ---- ##
  preprocess_recipe_obj <-
    sample_tbl |>
    recipe(fixed) |>
    update_role(smp_weights, new_role = "case_weight") |>
    update_role(domains, new_role = "id variable") |>
    step_string2factor(all_nominal_predictors()) |>
    step_other(all_nominal_predictors(), threshold = 0.01) |>
    prep() |>
    bake(new_data = NULL)

  ## ---- create folds ---- ##
  cvfold_obj <-
    rsample::group_vfold_cv(preprocess_recipe_obj,
                            group = cluster,
                            v = folds)

  ## ---- dynamic XGBoost model specification ---- ##
  # Helper function to set tune() or fixed
  # param_or_tune <- function(param) {
  #   if (param %in% names(tune_params)) tune() else if (!is.null(tune_params[[param]])) tune_params[[param]][1] else NULL
  # }

  xgboost_model <-
    boost_tree(mode = "regression",
               trees = define_param("trees", tune_params),
               mtry = define_param("mtry", tune_params),
               min_n = define_param("min_n", tune_params),
               tree_depth = define_param("tree_depth", tune_params),
               learn_rate = define_param("learn_rate", tune_params),
               loss_reduction = define_param("loss_reduction", tune_params),
               sample_size = define_param("sample_size", tune_params),
               stop_iter = define_param("stop_iter", tune_params)) %>%
    set_engine("xgboost",
               #objective = "reg:squarederror",
               alpha = define_param("alpha", tune_params),
               lambda = define_param("lambda", tune_params),
               colsample_bytree = define_param("colsample_bytree", tune_params),
               colsample_bylevel = define_param("colsample_bylevel", tune_params),
               colsample_bynode = define_param("colsample_bynode", tune_params),
               max_delta_step = define_param("max_delta_step", tune_params),
               weight = smp_weights,
               counts = FALSE)

  ### boost_tree uses tidy eval to enquo() its arguments, eval_enquo_xgbargs() has
  ### been created to evaluate boost_tree arguments by force and then enquo() the
  ### results of the evaluation.
  xgboost_model$args <- eval_enquo_xgbargs(xgboost_model$args) ## for the boost_tree args

  xgboost_model$eng_args <- eval_enquo_xgbargs(xgboost_model$eng_args) ## for set_engine args

  ## ---- dynamic parameter object ---- ##

  ## ---- create tuning grid ---- ##
  ### compute the appropriate grid size

  if (is.null(tune_size)){

    tune_size <- estimate_grid_size(tune_params = tune_params)

  }

  # Compute which parameters are actually set to tune()
  tuneable_params <- names(Filter(function(x) length(x) > 1, tune_params))

  # Only build param set for these
  xgb_param_list <- build_xgb_param_set(tune_params[tuneable_params])

  # Now create the grid
  xgboost_grid <- grid_space_filling(x = xgb_param_list, size = tune_size)

  ## ---- workflow ---- ##
  xgboost_wf <-
    workflows::workflow() |>
    workflows::add_model(xgboost_model) |>
    workflows::add_formula(fixed) |>
    workflows::add_case_weights(where(~smp_weights))

  ## quickly define the new metric i.e. domain specified rmse's computed with
  ## domain_rmse_vec() function
  domain_rmse <- yardstick::new_numeric_metric(
    domain_rmse_vec,
    direction = "minimize"
  )


  ## ---- tuning ---- ##
  if (!is.null(parallel_over)){

    # Detect number of available cores
    n_cores <- parallel::detectCores()

    # Set up parallel plan with multisession (safe for Windows)
    future::plan(multisession, workers = n_cores - 1)  # use all but one core


    ## tuning operation
    xgboost_tuned <-
      tune_grid(xgboost_wf,
                resamples = cvfold_obj,
                grid = xgboost_grid,
                metrics = metric_set(rmse, domain_rmse),
                control = control_grid(verbose = TRUE,
                                       save_pred = TRUE,
                                       parallel_over = parallel_over))

    # close the cores by enforcing 1 core
    future::plan(sequential)



  } else {

    xgboost_tuned <-
      tune_grid(xgboost_wf,
                resamples = cvfold_obj,
                grid = xgboost_grid,
                metrics = metric_set(rmse, domain_rmse),
                control = control_grid(verbose = TRUE))

  }

  ## ---- select best parameters ---- ##
  xgboost_best_params <-
    xgboost_tuned |>
    select_best(metric = "domain_rmse")


  return_obj <- list("xgb_tune_fullresults" = xgboost_tuned,
                     "xgb_best_params" = xgboost_best_params)

  class(return_obj) <- "xgb_tune_list"


  return(return_obj)

}



#' Determine a recommended tuning grid size
#' @param tune_params list of tuning parameters
#' @return An integer: recommended grid size
#'
estimate_grid_size <- function(tune_params) {
  # param_list should be a list of parameters (with NULLs if not tuned)
  p <- length(Filter(Negate(is.null), tune_params))  # count non-NULL params

  if (p <= 7) {
    size <- 2^p
  } else {
    size <- 17 * p
  }

  return(size)
}





build_xgb_param_set <- function(tune_params) {

  param_objs <- imap(tune_params, ~{
    param <- .y
    val   <- .x

    # skip NULLs
    if (is.null(val)) return(NULL)

    # fixed value (not tunable)
    if (length(val) == 1) return(val)

    # helper for generic numeric params
    new_quant_param_wrap <- function() {
      dials::new_quant_param(
        type = "double",
        range = val,
        inclusive = c(TRUE, TRUE)
      )
    }

      # tunable ranges
      switch(param,
             # built-in dials parameters
             trees             = trees(range = val),
             min_n             = min_n(range = val),
             tree_depth        = tree_depth(range = val),
             learn_rate        = learn_rate(range = val, trans = log10_trans()),
             loss_reduction    = loss_reduction(range = val),
             sample_size       = sample_prop(range = val),
             stop_iter         = stop_iter(range = val),
             alpha             = penalty(range = val, trans = log10_trans()),
             lambda            = penalty(range = val, trans = log10_trans()),

             # xgboost engine-specific params
             colsample_bytree  = new_quant_param_wrap(),
             colsample_bylevel = new_quant_param_wrap(),
             colsample_bynode  = new_quant_param_wrap(),
             max_delta_step    = new_quant_param_wrap(),

             # fallback (any other numeric)
             new_quant_param_wrap())

  })

  # drop NULLs and keep only "param" objects
  param_objs <- compact(param_objs)
  param_objs <- keep(param_objs, ~inherits(.x, "param"))

  # combine into a compact params object
  param_objs <- dials::parameters(x = param_objs)
  return(param_objs)
}


define_param <- function(param, tune_params) {

# If the parameter exists in the tune_params list AND has a vector (length>1) → tune()

  if (!is.null(tune_params[[param]]) && length(tune_params[[param]]) > 1) {
    return(tune())
  }

  # If the parameter exists in tune_params AND has a single value → fixed
  if (!is.null(tune_params[[param]]) && length(tune_params[[param]]) == 1) {
    return(tune_params[[param]][1])
  }

  # If the parameter is missing from tune_params → use default (NA)
  return(NULL)
}


eval_enquo_xgbargs <- function(args_list){

  xgbargs_list <-
  lapply(args_list,
         function(arg) {

           val <- eval_tidy(arg)              # evaluate the quosure expression
           env <- quo_get_env(arg)            # extract original environment
           result <- new_quosure(val, env = env)        # rewrap with preserved env

           return(result)
         })

  return(xgbargs_list)

}


domain_rmse_vec <- function(truth,
                            estimate,
                            na_rm = TRUE,
                            case_weights = NULL,
                            domains = NULL,
                            ...) {
  # if no weights provided, default to 1
  w <- case_weights %||% rep(1, length(truth))

  df_tbl <- tibble(truth = truth,
                   estimate = estimate,
                   domains = domains,
                   w = w)

  summary_tbl <- df_tbl %>%
    group_by(domains) %>%
    summarize(
      mean_truth = weighted.mean(truth, w = w, na.rm = na_rm),
      mean_est   = weighted.mean(estimate, w = w, na.rm = na_rm),
      sumw       = sum(w, na.rm = na_rm),
      .groups    = "drop"
    ) %>%
    summarize(rmse = sqrt(weighted.mean((mean_est - mean_truth)^2,
                                        w = sumw,
                                        na.rm = na_rm)))

  return(summary_tbl$rmse)
}














