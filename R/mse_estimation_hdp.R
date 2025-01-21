# Internal documentation -------------------------------------------------------

# MSE estimation - parametric bootstrap procedure

# Function parametric_bootstrap conducts the MSE estimation defined in function
# mse_estim_hdp (see below)
# The parametric boostrap approach can be find in Lahiri and Salvati (2023)

parametric_bootstrap_hdp <- function(framework,
                                     point_estim,
                                     fixed,
                                     transformation,
                                     interval = c(-1, 2),
                                     L,
                                     B,
                                     #boot_type,
                                     #parallel_mode,
                                     #cpus,
                                     benchmark,
                                     benchmark_type,
                                     benchmark_level) {


  message("\r", "Bootstrap started                                            ")

  res_s <- NULL
  fitted_s <- NULL

  start_time <- Sys.time()

  mses <- simplify2array(lapply(
    X = seq_len(B),
    FUN = mse_estim_wrapper_hdp,
    B = B,
    framework = framework,
    lambda = point_estim$optimal_lambda,
    shift = point_estim$shift_par,
    model_par = point_estim$model_par,
    gen_model = point_estim$gen_model,
    fixed = fixed,
    transformation = transformation,
    interval = interval,
    L = L,
    res_s = res_s,
    fitted_s = fitted_s,
    start_time = start_time,
    #boot_type = boot_type,
    benchmark = benchmark,
    benchmark_type = benchmark_type,
    benchmark_level = benchmark_level
  ))

  message("\r", "Bootstrap completed", "\n")
  if (.Platform$OS.type == "windows") {
    flush.console()
  }

  mses <- apply(mses, c(1, 2), mean)
  if(is.null(framework$aggregate_to_vec)){
    mses <- data.frame(Domain = unique(framework$pop_domains_vec), mses)
  }else{
    mses <- data.frame(Domain = unique(framework$aggregate_to_vec), mses)
  }

  return(mses)
}


# mse_estim_hdp (only internal) ----------------------------------------------------

# The mse_estim_hdp function defines all parameters and estimations which have to
# be replicated B times for the Parametric Bootstrap Approach.


mse_estim_hdp <- function(framework,
                          lambda,
                          shift,
                          model_par,
                          gen_model,
                          res_s,
                          fitted_s,
                          fixed,
                          transformation,
                          interval,
                          L,
                          #boot_type,
                          benchmark,
                          benchmark_type,
                          benchmark_level) {
  # The function superpopulation returns an income vector and a temporary
  # variable that passes the random effect to generating bootstrap populations
  # in bootstrap_par_hdp.
  superpop <- superpopulation_hdp(
    framework = framework,
    model_par = model_par,
    gen_model = gen_model,
    lambda = lambda,
    shift = shift,
    transformation = transformation,
    fixed = fixed
  )

  pop_income_vector <- superpop$pop_income_vector

  if (inherits(framework$threshold, "function")) {
    framework$threshold <-
      framework$threshold(y = pop_income_vector)
  }

  if(!is.null(framework$aggregate_to_vec)) {
    N_dom_pop_tmp <- framework$N_dom_pop_agg
    pop_domains_vec_tmp <- framework$aggregate_to_vec
  } else {
    N_dom_pop_tmp <- framework$N_dom_pop
    pop_domains_vec_tmp <- framework$pop_domains_vec
  }

  if(!is.null(framework$pop_weights)) {
    pop_weights_vec <- framework$pop_data[[framework$pop_weights]]
  }else{
    pop_weights_vec <- rep(1, nrow(framework$pop_data))
  }

  # True indicator values
  true_indicators <- matrix(
    nrow = N_dom_pop_tmp,
    data = unlist(lapply(framework$indicator_list,
                         function(f, threshold) {
                           matrix(
                             nrow = N_dom_pop_tmp,
                             data =
                               unlist(mapply(
                                 y = split(pop_income_vector, pop_domains_vec_tmp),
                                 pop_weights = split(pop_weights_vec, pop_domains_vec_tmp),
                                 f,
                                 threshold = framework$threshold
                               )),
                             byrow = TRUE
                           )
                         },
                         threshold = framework$threshold
    ))
  )

  colnames(true_indicators) <- framework$indicator_names

  if (!is.null(benchmark)) {
    if (is.character(benchmark)) {
      add_bench <- data.frame(true_indicators[, benchmark])
      if (!is.null(dim(add_bench))) {
        colnames(add_bench) <- c(paste0(benchmark,"_bench"))
      } else {
        names(add_bench) <- c(paste0(benchmark,"_bench"))
      }
    } else {
      if (is.numeric(benchmark)) {
        add_bench <- true_indicators[, names(benchmark)]
        if (!is.null(dim(add_bench))) {
          colnames(add_bench) <- c(paste0(names(benchmark),"_bench"))
        } else {
          names(add_bench) <- c(paste0(names(benchmark),"_bench"))
        }
      } else {
        add_bench <- true_indicators[, names(benchmark)[-1]]
        if (!is.null(dim(add_bench))) {
          colnames(add_bench) <- c(paste0(names(benchmark)[-1],"_bench"))
        } else {
          names(add_bench) <- c(paste0(names(benchmark)[-1],"_bench"))
        }
      }
    }
    true_indicators <- as.matrix(cbind(true_indicators, add_bench))
  }

  # The function bootstrap_par_hdp returns a sample that can be given into the
  # point estimation to get predictors of the indicators that can be compared
  # to the "truth".

  bootstrap_sample <- bootstrap_par_hdp(
    fixed = fixed,
    transformation = transformation,
    framework = framework,
    model_par = model_par,
    lambda = lambda,
    shift = shift,
    vu_tmp = superpop$vu_tmp
  )


  framework$smp_data <- bootstrap_sample

  # Prediction of indicators with bootstap sample.
  bootstrap_point_estim <- as.matrix(point_estim_hdp(
    fixed = fixed,
    transformation =
      transformation,
    interval = interval,
    L = L,
    framework = framework
  )[[1]][, -1])

  # benchmark
  if (!is.null(benchmark)) {
    if (is.null(benchmark_level)) {
      bootstrap_point_estim <- benchmark_ebp_national(
        point_estim = bootstrap_point_estim,
        framework = framework,
        fixed = fixed,
        benchmark = benchmark,
        benchmark_type = benchmark_type)
    } else {
      bootstrap_point_estim <- benchmark_ebp_level(
        point_estim = bootstrap_point_estim,
        framework = framework,
        fixed = fixed,
        benchmark = benchmark,
        benchmark_type = benchmark_type,
        benchmark_level = benchmark_level)
    }
  }
  return((bootstrap_point_estim - true_indicators)^2)
} # End mse_estim


# Superpopulation function -----------------------------------------------------

# The model parameter from the nested error linear regression model are
# used to contruct a superpopulation model.
superpopulation_hdp <- function(framework, model_par, gen_model, lambda, shift,
                                transformation, fixed) {
  # superpopulation individual errors
  eps <- rep(0, framework$N_pop)
  eps[framework$obs_dom] <- rnorm(sum(framework$n_pop[framework$dist_obs_dom]),
                                  0,
                                  rep(sqrt(model_par$sigmae2est), framework$n_pop[framework$dist_obs_dom]))
  if(framework$N_dom_unobs>0){
  eps[!framework$obs_dom] <- rnorm(sum(framework$n_pop[!framework$dist_obs_dom]),
                                       0,
                                       rep(sqrt(model_par$sigma2e.out), framework$n_pop[!framework$dist_obs_dom]))
    }
  # superpopulation random effect
  vu_tmp <- rnorm(framework$N_dom_pop, 0, sqrt(model_par$sigmau2est))
  vu_pop <- rep(vu_tmp, framework$n_pop)

  #  superpopulation income vector
  Y_pop_b <- gen_model$mu_fixed + eps + vu_pop

  Y_pop_b <- back_transformation(
    y = Y_pop_b,
    transformation = transformation,
    lambda = lambda,
    shift = shift,
    framework = framework,
    fixed = fixed
  )
  Y_pop_b[!is.finite(Y_pop_b)] <- 0

  return(list(pop_income_vector = Y_pop_b, vu_tmp = vu_tmp))
}

# Bootstrap function -----------------------------------------------------------

bootstrap_par_hdp <- function(fixed, transformation, framework, model_par, lambda,
                              shift, vu_tmp) {
  # Bootstrap sample individual error term
  eps <- rnorm(framework$N_smp, 0, rep(sqrt(model_par$sigmae2est), framework$n_smp))

  # Bootstrap sample random effect
  vu_smp <- rep(vu_tmp[framework$dist_obs_dom], framework$n_smp)
  # Extraction of design matrix
  X_smp <- model.matrix(fixed, framework$smp_data)
  # Constant part of income vector for bootstrap sample
  betas_in <- apply(model_par$betas, 1, rep, framework$n_smp)
  mu_smp <- diag(X_smp %*% t(betas_in))
  # Transformed bootstrap income vector
  Y_smp_b <- mu_smp + eps + vu_smp
  # Back transformation of bootstrap income vector
  Y_smp_b <- back_transformation(
    y = Y_smp_b,
    transformation = transformation,
    lambda = lambda,
    shift = shift,
    framework = framework,
    fixed = fixed
  )
  Y_smp_b[!is.finite(Y_smp_b)] <- 0

  # Inclusion of bootstrap income vector into sample data
  bootstrap_smp <- framework$smp_data
  bootstrap_smp[paste(fixed[2])] <- Y_smp_b

  return(bootstrap_sample = bootstrap_smp)
}


# progress for mse_estim_hdp (only internal) ----------
mse_estim_wrapper_hdp <- function(i,
                                  B,
                                  framework,
                                  lambda,
                                  shift,
                                  model_par,
                                  gen_model,
                                  fixed,
                                  transformation,
                                  interval,
                                  L,
                                  res_s,
                                  fitted_s,
                                  start_time,
                                  # boot_type,
                                  seedvec,
                                  benchmark,
                                  benchmark_type,
                                  benchmark_level) {

  tmp <- mse_estim_hdp(
    framework = framework,
    lambda = lambda,
    shift = shift,
    model_par = model_par,
    gen_model = gen_model,
    res_s = res_s,
    fitted_s = fitted_s,
    fixed = fixed,
    transformation = transformation,
    interval = interval,
    L = L,
    # boot_type = boot_type,
    benchmark = benchmark,
    benchmark_type = benchmark_type,
    benchmark_level = benchmark_level
  )

  if (i %% 10 == 0) {
    if (i != B) {
      delta <- difftime(Sys.time(), start_time, units = "secs")
      remaining <- (delta / i) * (B - i)
      remaining <- unclass(remaining)
      remaining <- sprintf(
        "%02d:%02d:%02d:%02d",
        remaining %/% 86400, # days
        remaining %% 86400 %/% 3600, # hours
        remaining %% 3600 %/% 60, # minutes
        remaining %% 60 %/% 1
      ) # seconds)

      message("\r", i, " of ", B, " Bootstrap iterations completed \t
              Approximately ", remaining, " remaining \n")
      if (.Platform$OS.type == "windows") flush.console()
    }
  }
  return(tmp)
}
