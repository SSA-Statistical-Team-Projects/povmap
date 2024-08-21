#' Empirical Best Prediction for Disaggregated Indicators 
#' under A Nested Error Regression Model with High Dimensional Parameter 
#' @param fixed a two-sided linear formula object describing the
#' fixed-effects part of the nested error linear regression model with the
#' dependent variable on the left of a ~ operator and the explanatory
#' variables on the right, separated by + operators. The argument corresponds
#' to the argument \code{fixed} in function \code{\link[nlme]{lme}}.
#' @param pop_data a data frame that needs to comprise the variables
#' named on the right of the ~ operator in \code{fixed}, i.e. the explanatory
#' variables, and \code{pop_domains}.
#' @param pop_domains a character string containing the name of a variable that
#' indicates domains in the population data. The variable can be numeric or
#' a factor but needs to be of the same class as the variable named in
#' \code{smp_domains}.
#' @param smp_data a data frame that needs to comprise all variables named in
#' \code{fixed} and \code{smp_domains}.
#' @param smp_domains a character string containing the name of a variable
#' that indicates domains in the sample data. The variable can be numeric or a
#' factor but needs to be of the same class as the variable named in
#' \code{pop_domains}.
#' @param threshold a number defining a threshold. Alternatively, a threshold
#' may be defined as a \code{function} of \code{y} returning a numeric value.
#' Such a function will be evaluated once for the point estimation and in each
#' iteration of the parametric bootstrap. A threshold is needed for calculation
#' e.g. of head count ratios and poverty gaps. The  argument defaults to
#' \code{NULL}. In this case, the threshold is set to 60\% of the median of the
#' variable that is selected as dependent variable similary to the
#' at-risk-of-poverty rate used in the EU (see also
#' \cite{Social Protection  Committee 2001}). However, any desired threshold can
#' be chosen.
#' @param transformation a character string. Five different transformation
#' types for the dependent variable can be chosen (i) no transformation ("no");
#' (ii) log transformation ("log"); (iii) Box-Cox transformation ("box.cox");
#' (iv) Dual transformation ("dual"); (v) Log-Shift transformation
#' ("log.shift"); (Vi) rank-order transformation ("ordernorm"). Defaults to
#' \code{"box.cox"}.
#' @param interval a string equal to 'default' or a numeric vector containing a
#' lower and upper limit determining an interval for the estimation of the
#' optimal parameter. The interval is passed to function
#' \code{\link[stats]{optimize}} for the optimization. Defaults to 'default'
#' which equals c(-1,2) for Box-Cox, c(0,2) for Dual and an interval based on
#' the range of y for Log-Shift transformation. If the convergence fails, it is
#' often advisable to choose a smaller more suitable interval. For right skewed
#' distributions, the negative values may be excluded, also values larger than
#' 1 are seldom observed.
#' @param L a number determining the number of Monte-Carlo simulations that
#' must be at least 1. Defaults to 50. For practical applications, values
#' larger than 200 are recommended (see also
#' \cite{Molina, I. and Rao, J.N.K. (2010)}).
#' @param MSE if \code{TRUE}, MSE estimates using a parametric bootstrap
#' approach are calculated (see also \cite{Gonzalez-Manteiga et al. (2008)}).
#' Defaults to \code{FALSE}.
#' @param B a number determining the number of bootstrap populations in the
#' parametric bootstrap approach (see also
#' \cite{Gonzalez-Manteiga et al. (2008)}) used in the MSE estimation. The
#' number must be greater than 1. Defaults to 50. For practical applications,
#' values larger than 200 are recommended (see also
#' \cite{Molina, I. and Rao, J.N.K. (2010)}).
#' @param na.rm if \code{TRUE}, observations with \code{NA} values are deleted
#' from the population and sample data. For the EBP procedure complete
#' observations are required. Defaults to \code{FALSE}.
#' @param weights a character string containing the name of a variable that
#' indicates weights in the sample data. If a character string is provided
#' a weighted version of the ebp will be used. The variable has to be numeric.
#' Defaults to \code{NULL}.
#' @param pop_weights a character string containing the name of a variable that
#' indicates population weights in the populatation data. If a character string
#' is provided weighted indicators are estimated using population weights.
#' The variable has to be numeric. Defaults to \code{NULL}.
#' @param aggregate_to a character string containing the name of a variable from
#' population data that indicates the target domain level for which the
#' results are to be displayed. The variable can be numeric or a factor.
#' Defaults to \code{NULL}.
#' @param benchmark The input depends on the type of benchmarking to be
#' performed.
#' (i) Benchmarking with a fixed value:
#' (a) with one value for each indicator: a named vector containing the numeric
#' benchmark value(s). The names of the vector matchs to the chosen indicators.
#' Benchmarking is available for \code{"Mean"} and \code{"Head_Count"}.
#' (b) with values for the sub-level specified in the argument
#' \code{benchmark_level}: a data.frame composed of a variable of class
#' character containing the domain names at which the benchmarkaing is
#' performed and variable(s) with benchmark value(s) of class numeric.
#' Benchmarking is supplied for the Mean and the Head_Count ratio. Therefore,
#' the names of the data.frame must match for the first variable the
#' benchmark_level and for the other(s) to Mean and Head_Count.
#' (ii) Benchmarking with the survey data: a vector containing the names of the
#' chosen indicators. In this case, survey weights (\code{weights}) are needed.
#' Benchmarking is available for \code{"Mean"} and \code{"Head_Count"}.
#' @param benchmark_type a character indicating the type of benchmarking. Types
#' that can be chosen (i) Raking ("\code{raking}"), (ii) Ratio adjustment
#' ("\code{ratio}"), and for head count, ratio adjustment of the complement
#' ("\code{ratio_complement}". Defaults to "\code{ratio}"
#' @param benchmark_weights the name of variable containing benchmark weights.
#' This is only possible for internal benchmarking and enable users to benchmark
#' with weights differing from the survey weights (Default for weighting for
#' internal benchmarking).
#' @param Q the quantiles for fitting the M-quantile regression.
#' @param method the method is used to obtain area-specific tuning parameter. 
#' Two methods can be chosen
#' (i) "LBP": linear best predictor (ii) "Naive": traditional method.
#' @param k_sigma_u tuning constant used for Huber influence function for 
#' the estimation of the sigma2u.
#' @param k_sigma_e tuning constant used for Huber influence function for 
#' the estimation of the sigma2ei.
#' @param k_b=1.345 tuning constant used for Huber influence function for 
#' the estimation of the area-specific coefficients.
ebp_hdp <- function(fixed,
                    pop_data,
                    pop_domains,
                    smp_data,
                    smp_domains,
                    L = 50,
                    threshold = NULL,
                    transformation = "box.cox",
                    interval = "default",
                    MSE = FALSE,
                    B = 50,
                    na.rm = FALSE,
                    weights = NULL,
                    pop_weights = NULL,
                    aggregate_to = NULL,
                    custom_indicator = NULL,
                    benchmark = NULL,
                    benchmark_type = "ratio",
                    benchmark_weights = NULL,
                    Q = sort(c(seq(0.006,0.99,0.045),0.5,0.994,0.01,0.02,0.96,0.98)),
                    method=c("LBP","Naive"),
                    tol=1e-06,
                    maxit=100,
                    k_b=1.345,
                    k_sigma_u=2.5,
                    k_sigma_e=2.5
                    ){
  # Check the arguments --------------------------------------------------------
  ebphdp_check1(
    fixed = fixed, pop_data = pop_data, pop_domains = pop_domains,
    smp_data = smp_data, smp_domains = smp_domains, L = L
  )
  
  ebphdp_check2(
    threshold = threshold, 
    transformation = transformation,
    interval = interval,
    MSE = MSE,
    B = B,
    L = L,
    custom_indicator = custom_indicator,
    na.rm = na.rm,
    weights = weights,
    pop_weights = pop_weights,
    benchmark = benchmark,
    benchmark_type = benchmark_type,
    benchmark_level = benchmark_level,
    benchmark_weights = benchmark_weights
  )
  
  # Save function call ---------------------------------------------------------
  
  call <- match.call()
  if (inherits(call$fixed, "name")) {
    call$fixed <- fixed
  }
  
    # Data manipulation and notational framework ---------------------------------

  if (is.null(benchmark_weights) & !is.null(weights)) {
    benchmark_weights <- "benchmark_weights"
    smp_data$benchmark_weights <- smp_data[[weights]]
  }

  # The function framework_HDP can be found in script framework_HDP.R
  framework <- framework_hdp(
    fixed = fixed,
    pop_data = pop_data,
    pop_domains = pop_domains,
    smp_data = smp_data,
    smp_domains = smp_domains,
    threshold = threshold,
    na.rm = na.rm,
    aggregate_to = aggregate_to,
    weights = weights,
    pop_weights = pop_weights,
    custom_indicator = custom_indicator,
    Q = Q,
    method = method,
    tol = tol, 
    maxit = maxit,
    k_b = k_b,
    #k_sigma = k_sigma,
    k_sigma_u = k_sigma_u,
    k_sigma_e = k_sigma_e
  )
  
  # Point Estimation -----------------------------------------------------------
  point_estim <- point_estim_hdp(
    framework = framework,
    fixed = fixed,
    transformation = transformation,
    interval = interval,
    L = L,
    keep_data = TRUE,
    Ydump = Ydump
  )
  
  # benchmarking
  if (!is.null(benchmark)) {
    if (is.null(benchmark_level)) {
      point_estim$ind <- benchmark_ebp_national(
        point_estim = point_estim,
        framework = framework,
        fixed = fixed,
        benchmark = benchmark,
        benchmark_type = benchmark_type)
    } else {
      point_estim$ind <- benchmark_ebp_level(
        point_estim = point_estim,
        framework = framework,
        fixed = fixed,
        benchmark = benchmark,
        benchmark_type = benchmark_type,
        benchmark_level = benchmark_level)
    }
  
    if (any(names(point_estim$ind) %in% c("Mean_bench"))) {
      if (any(is.na(point_estim$ind$Mean_bench))) {
        message(strwrap(prefix = " ", initial = "",
                        "Benchmark point estimates for
                          Mean contain missing values. Please check source data"))
      }
    }
    
    if (any(names(point_estim$ind) %in% c("Head_Count_bench"))) {
      if (any(is.na(point_estim$ind$Head_Count_bench))) {
        message(strwrap(prefix = " ", initial = "",
                        "Benchmark point estimates for
                          Head_Count contain missing values. Please check source data"))
      }
      else if(!all(point_estim$ind$Head_Count_bench >= 0 &
                   point_estim$ind$Head_Count_bench <= 1)){
        message(strwrap(prefix = " ", initial = "",
                        "Please note that benchmark point estimates for
                          Head_Count are outside the expected range [0,1]."))
      }
    }
  }
  
  # MSE Estimation -------------------------------------------------------------
  if (MSE == TRUE) {
    # The function parametric_bootstrap can be found in script mse_estimation_hdp.R
    mse_estimates <- parametric_bootstrap_hdp(
      framework = framework,
      point_estim = point_estim,
      fixed = fixed,
      transformation = transformation,
      interval = interval,
      L = L,
      B = B,
      # boot_type = boot_type,
      # parallel_mode = parallel_mode,
      # cpus = cpus,
      benchmark = benchmark,
      benchmark_type = benchmark_type,
      benchmark_level = benchmark_level
    )
    
    
    hdp_out <- list(
      ind = point_estim$ind,
      MSE = mse_estimates,
      transform_param = point_estim[c(
        "optimal_lambda",
        "shift_par"
      )],
      model = point_estim$model,
      model_par = point_estim$model_par,
      framework = framework[c(
        "N_dom_unobs",
        "N_dom_smp",
        "N_smp",
        "N_pop",
        "smp_domains",
        "smp_data",
        "smp_domains_vec",
        "pop_domains_vec"
      )],
      transformation = transformation,
      method = "reml",
      fixed = fixed,
      call = call,
      successful_bootstraps = NULL
    )
  } 
  else {
  hdp_out <- list(
    ind = point_estim$ind,
    MSE = NULL,
    transform_param = point_estim[c(
      "optimal_lambda",
      "shift_par"
    )],
    model = point_estim$model,
    model_par = point_estim$model_par,
    framework = framework[c(
      "N_dom_unobs",
      "N_dom_smp",
      "N_smp",
      "N_pop",
      "smp_domains",
      "smp_data",
      "smp_domains_vec",
      "pop_domains_vec"
    )],
    transformation = transformation,
    method = "reml",
    fixed = fixed,
    call = call,
    successful_bootstraps = NULL
  )
  }
  
  class(hdp_out) <- c("ebphdp", "emdi")
  return(hdp_out)
}
