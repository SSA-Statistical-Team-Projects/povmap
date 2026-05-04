
#' MSE Estimation for Mixed Effect Gradient Boosting
#'
#' This function calculates the Mean Squared Error (MSE) Estimates for a
#' Mixed Effect Gradient Boosting model while using bootstrap samples.
#'
#' @param Y numeric, target variable
#' @param X data.frame, Covariates
#' @param dom_name factor, name of domain variable
#' @param smp_data data.frame, sample data
#' @param model list, trained MEGB model
#' @param error_sd numeric, standard deviation of error
#' @param pop_data data.frame, population data
#' @param B integer, number of bootstrap samples for MSE
#' @param initial_random_effects numeric, initial value for Random Effects
#' @param ErrorTolerance numeric, tolerance for EM-algorithm
#' @param MaxIterations integer, maximal number of iterations for the EM-algorithm
#' @param cov_names character, names of covariates
#' @param gradient_params list, optional list with Gradient Boosting parameters
#' @param formula_random_effects formula, formula for Random Effects
#' @param bootstrap_cores integer, number of cores to use for bootstrap
#' @param seed integer, choose seed
#' @param ... additional arguments
#' @import future
#' @importFrom furrr future_map
#' @importFrom purrr map
#' @importFrom xgboost xgboost xgb.DMatrix xgb.importance
#' @importFrom dplyr group_by summarise
#' @return list of results of MSE estimates


# mse_megb -------------------------------------------------------------------------


mse_megb <- function(Y,                           # target variable
                     X,                           # covariates
                     dom_name,                    # name of domain variable
                     smp_data,                    # sample data
                     model,                       # trained Mixed Effect GB model
                     error_sd,                    # standard deviation of Errors
                     pop_data,                    # population data
                     B,                     # number of bootstrap samples
                     initial_random_effects,    # initial value for Random Effects
                     ErrorTolerance,     # error tolerance for EM-algorithm
                     MaxIterations,          # max number of EM iterations
                     cov_names,                   # names of covariats
                     gradient_params = list(),    # Gradient Boosting parameters
                     formula_random_effects, # formula for Random Effects
                     bootstrap_cores,    # number of bootstrap cores
                     seed,
                     gbm_engine,
                     unit_pred_smp,
                     unit_preds,
                     corrected_bt = NULL,  # Jensen-corrected back-transform function
                     benchmark_fn = NULL,  # function(est_orig, boot_smp_data) -> bench
                     ...) {

  # extract domains and check if part of sample data
  domains <- t(unique(pop_data[dom_name]))
  in_samp <- domains %in% t(unique(smp_data[dom_name]))

  # number of observations per domain
  n_i <- as.numeric(table(pop_data[[dom_name]]))

  # prepare data for sampling
  ran_obj <- ran_comp(
    Y = Y,
    smp_data = smp_data,
    unit_pred_smp = unit_pred_smp,
    error_sd = error_sd,
    dom_name = dom_name,
    cov_names = cov_names,
    model = model
  )
  ran_effs <- ran_obj$ran_effs
  gb_res <- ran_obj$gb_res
  smp_data <- ran_obj$smp_data

  # generate matrix of predicted population data
  pred_mat <- matrix(unit_preds$unit_preds, nrow = length(unit_preds$unit_preds), ncol = B)

  # sampling-function for block-sampling of errors
  block_sample_e <- function(x) {
    block_sample(
      domains = domains,
      in_samp = in_samp,
      smp_data = smp_data,
      dom_name = dom_name,
      pop_data = pop_data,
      gb_res = gb_res
    )
  }

  # sampling-function for Random Effects
  sample_ui <- function(x) {
    rep(
      sample(ran_effs, size = length(n_i), replace = TRUE),
      n_i
    )
  }

  # initialise matrix for residual errors
  e_ij <- matrix(NA, nrow = length(unit_preds$unit_preds), ncol = B)
  e_ij <- apply(e_ij, 2, block_sample_e)

  # calculate random effects for bootstrap
  u_i <- apply(pred_mat, 2, sample_ui)

  # remove unnecessary column
  smp_data$gb_res <- NULL

  # combine predictions, Random Effects and error to form y_star
  y_star <- pred_mat + u_i + e_ij

  # group indices of observations by domain
  indi_agg <- rep(1:length(n_i), n_i)

  # function to aggregate for domains
  my_agg <- function(x) {
    tapply(x, indi_agg, mean)
  }

  # Calculate tau* (aggregated bootstrap predictions)
  tau_star <- apply(y_star, MARGIN = 2, my_agg)

  # create bootstrap samples
  boots_sample <- vector(mode = "list", length = B)
  for (i in 1:B) {
    pop_data$y_star <- y_star[, i]
    boots_sample[[i]] <- sample_select(pop_data, smp = smp_data, dom_name = dom_name)
  }

  # predict function for bootstrap samples
  my_estim_f <- function(x) {
    model_boot <- em_gb_lmm(
              Y = x$y_star,
              X = x[, colnames(X)],
              dom_name = dom_name,
              data = x,
              pop_data = pop_data,
              seed = 1,
              gradient_params = gradient_params,
              formula_random_effects = formula_random_effects,
              initial_random_effects = 0, # initial values for domain-effects
              max_iterations = 10,        # max EM-iterations
              error_tolerance = 1e-04,
              cov_names = cov_names,
              gbm_engine = gbm_engine,
              ...
    )

    unit_level_predictions <- gbm_predict(
      model = model_boot,
      smp_data = smp_data,
      pop_data = pop_data,
      Y = Y,
      dom_name = dom_name,
      gbm_engine = gbm_engine,
      cov_names = cov_names
    )
    
    mean_preds <- unit_level_predictions$unit_pred_pop %>%
      group_by(dom_name) %>%
      dplyr::summarise(Mean = mean(unit_preds)) %>%
      as.data.frame()

    mean_boot_t <- mean_preds[, "Mean"]   # transformed-space domain means

    # Back-transform and benchmark if functions were supplied
    mean_boot_orig  <- if (!is.null(corrected_bt)) corrected_bt(mean_boot_t)
                       else mean_boot_t
    mean_boot_bench <- if (!is.null(benchmark_fn))
                         benchmark_fn(mean_boot_orig, x)
                       else NULL

    list(
      Mean_boot       = mean_boot_t,
      Mean_boot_bench = mean_boot_bench,
      error_sd_boot   = model_boot$error_sd,
      ran_eff_sd_boot = model_boot$ran_eff_sd
    )

  }

  # calculate tau_b (bootstrap estimates)
  if (is.null(bootstrap_cores) || bootstrap_cores <= 1) {
  boots_models <- purrr::map(boots_sample, my_estim_f, .progress = TRUE)
  }
  else{
    # recognise used system
    os_type <- Sys.info()[["sysname"]]

    if (os_type == "Windows") {
      future::plan(future::multisession, workers = bootstrap_cores)
    } else {
      future::plan(future::multicore, workers = bootstrap_cores)
    }
    boots_models <- furrr::future_map(
      boots_sample,
      my_estim_f,
      .options = furrr::furrr_options(seed = seed),
      .progress = TRUE
    )
    future::plan(future::sequential)
  }

  tau_b <- sapply(boots_models, getElement, "Mean_boot") |> unlist()
  # calculate MSE estimates (in transformed space; delta-method applied in wrapper)
  MSE_estimates <- rowMeans((tau_star - tau_b)^2)
  MSE_estimates <- data.frame(
    unique(pop_data[dom_name]),
    Mean = MSE_estimates
  )
  rownames(MSE_estimates) <- NULL

  # Benchmarked MSE in original space (no delta-method needed)
  MSE_bench_estimates <- NULL
  if (!is.null(corrected_bt) && !is.null(benchmark_fn)) {
    tau_star_orig <- apply(tau_star, 2, corrected_bt)   # ndomains x B matrix
    tau_b_bench   <- sapply(boots_models, getElement, "Mean_boot_bench")

    if (!is.null(tau_b_bench) && !all(sapply(tau_b_bench, is.null))) {
      mse_bench_vec <- rowMeans((tau_star_orig - tau_b_bench)^2, na.rm = TRUE)
      MSE_bench_estimates <- data.frame(
        unique(pop_data[dom_name]),
        Mean = mse_bench_vec
      )
      rownames(MSE_bench_estimates) <- NULL
    }
  }

  boot_error_sd <- sapply(boots_models, getElement, "error_sd_boot") |> unlist()
  boot_ran_eff_sd_boot <- sapply(boots_models, getElement, "ran_eff_sd_boot") |> unlist()
  # list of results (output)
  list(
    call = call,
    MSE_estimates       = MSE_estimates,
    MSE_bench_estimates = MSE_bench_estimates,
    boot_error_sd = boot_error_sd,
    boot_ran_eff_sd_boot = boot_ran_eff_sd_boot,
    error_sd_input = error_sd
  )
}


