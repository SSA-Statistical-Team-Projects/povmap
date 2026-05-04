# Internal: parametric bootstrap MSE for megb_em
#' @importFrom purrr map
#' @importFrom furrr future_map furrr_options
#' @importFrom future plan multisession multicore sequential
#' @importFrom dplyr group_by summarise

mse_megb <- function(Y, X, dom_name, smp_data, model, error_sd, pop_data,
                     B, initial_random_effects, ErrorTolerance, MaxIterations,
                     cov_names, gradient_params = list(), formula_random_effects,
                     bootstrap_cores, seed, gbm_engine,
                     unit_pred_smp, unit_preds,
                     corrected_bt = NULL,
                     benchmark_fn = NULL,
                     ...) {

  domains  <- t(unique(pop_data[dom_name]))
  in_samp  <- domains %in% t(unique(smp_data[dom_name]))
  n_i      <- as.numeric(table(pop_data[[dom_name]]))

  ran_obj  <- ran_comp(
    Y = Y, smp_data = smp_data, unit_pred_smp = unit_pred_smp,
    error_sd = error_sd, dom_name = dom_name,
    cov_names = cov_names, model = model
  )
  ran_effs <- ran_obj$ran_effs
  gb_res   <- ran_obj$gb_res
  smp_data <- ran_obj$smp_data

  pred_mat <- matrix(unit_preds$unit_preds,
                     nrow = length(unit_preds$unit_preds), ncol = B)

  block_sample_e <- function(x) {
    block_sample(domains = domains, in_samp = in_samp, smp_data = smp_data,
                 dom_name = dom_name, pop_data = pop_data, gb_res = gb_res)
  }

  sample_ui <- function(x) {
    rep(sample(ran_effs, size = length(n_i), replace = TRUE), n_i)
  }

  e_ij    <- matrix(NA, nrow = length(unit_preds$unit_preds), ncol = B)
  e_ij    <- apply(e_ij, 2, block_sample_e)
  u_i     <- apply(pred_mat, 2, sample_ui)
  smp_data$gb_res <- NULL

  y_star  <- pred_mat + u_i + e_ij
  indi_agg <- rep(1:length(n_i), n_i)
  my_agg   <- function(x) tapply(x, indi_agg, mean)
  tau_star <- apply(y_star, MARGIN = 2, my_agg)

  boots_sample <- vector(mode = "list", length = B)
  for (i in 1:B) {
    pop_data$y_star   <- y_star[, i]
    boots_sample[[i]] <- sample_select(pop_data, smp = smp_data, dom_name = dom_name)
  }

  my_estim_f <- function(x) {
    model_boot <- em_gb_lmm(
      Y                      = x$y_star,
      X                      = x[, colnames(X)],
      dom_name               = dom_name,
      data                   = x,
      pop_data               = pop_data,
      seed                   = 1,
      gradient_params        = gradient_params,
      formula_random_effects = formula_random_effects,
      initial_random_effects = 0,
      max_iterations         = 10,
      error_tolerance        = 1e-04,
      cov_names              = cov_names,
      gbm_engine             = gbm_engine,
      ...
    )

    unit_level_predictions <- gbm_predict(
      model      = model_boot,
      smp_data   = smp_data,
      pop_data   = pop_data,
      Y          = Y,
      dom_name   = dom_name,
      gbm_engine = gbm_engine,
      cov_names  = cov_names
    )

    mean_preds <- unit_level_predictions$unit_pred_pop |>
      dplyr::group_by(dom_name) |>
      dplyr::summarise(Mean = mean(unit_preds)) |>
      as.data.frame()

    mean_boot_t <- mean_preds[, "Mean"]

    mean_boot_orig  <- if (!is.null(corrected_bt)) corrected_bt(mean_boot_t)
                       else mean_boot_t
    mean_boot_bench <- if (!is.null(benchmark_fn)) benchmark_fn(mean_boot_orig, x)
                       else NULL

    list(
      Mean_boot       = mean_boot_t,
      Mean_boot_bench = mean_boot_bench,
      error_sd_boot   = model_boot$error_sd,
      ran_eff_sd_boot = model_boot$ran_eff_sd
    )
  }

  if (is.null(bootstrap_cores) || bootstrap_cores <= 1) {
    boots_models <- vector("list", B)
    for (i in seq_len(B)) {
      message("Bootstrap iteration ", i, " of ", B)
      boots_models[[i]] <- my_estim_f(boots_sample[[i]])
    }
  } else {
    os_type <- Sys.info()[["sysname"]]
    if (os_type == "Windows") {
      future::plan(future::multisession, workers = bootstrap_cores)
    } else {
      future::plan(future::multicore, workers = bootstrap_cores)
    }
    old_max <- getOption("future.globals.maxSize")
    options(future.globals.maxSize = 4 * 1024^3)  # 4 GiB: pop_data can be large
    tryCatch(
      boots_models <- furrr::future_map(
        boots_sample, my_estim_f,
        .options  = furrr::furrr_options(seed = seed),
        .progress = TRUE
      ),
      finally = {
        options(future.globals.maxSize = old_max)
        future::plan(future::sequential)
      }
    )
  }

  tau_b         <- sapply(boots_models, getElement, "Mean_boot") |> unlist()
  MSE_estimates <- rowMeans((tau_star - tau_b)^2)
  MSE_estimates <- data.frame(unique(pop_data[dom_name]), Mean = MSE_estimates)
  rownames(MSE_estimates) <- NULL

  MSE_bench_estimates <- NULL
  if (!is.null(corrected_bt) && !is.null(benchmark_fn)) {
    tau_star_orig <- apply(tau_star, 2, corrected_bt)
    tau_b_bench   <- sapply(boots_models, getElement, "Mean_boot_bench")
    if (!is.null(tau_b_bench) && !all(sapply(tau_b_bench, is.null))) {
      mse_bench_vec <- rowMeans((tau_star_orig - tau_b_bench)^2, na.rm = TRUE)
      MSE_bench_estimates <- data.frame(unique(pop_data[dom_name]), Mean = mse_bench_vec)
      rownames(MSE_bench_estimates) <- NULL
    }
  }

  boot_error_sd        <- sapply(boots_models, getElement, "error_sd_boot")   |> unlist()
  boot_ran_eff_sd_boot <- sapply(boots_models, getElement, "ran_eff_sd_boot") |> unlist()

  list(
    call                 = match.call(),
    MSE_estimates        = MSE_estimates,
    MSE_bench_estimates  = MSE_bench_estimates,
    boot_error_sd        = boot_error_sd,
    boot_ran_eff_sd_boot = boot_ran_eff_sd_boot,
    error_sd_input       = error_sd
  )
}
