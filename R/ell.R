#' Prediction for Disaggregated Indicators Using ELL 
#'
#' Function \code{ebp} estimates indicators using the
#' approach by \cite{Elbers, Lanjouw, and Lanjouw (2010)}. Point predictions of indicators
#' are obtained by Monte-Carlo approximations based on an unconditional random effect 
#' model. Additionally, variance is estimated treating parameters as uncertain. 
#' As in EBP, five different transformation types for the dependent variable can be chosen.
#' This approach can be extended to data under informative sampling using
#' weights and is based on \cite{Guadarrama et al. (2018)}.
#'
#' @param fixed a two-sided linear formula object describing the
#' linear regression model with the dependent variable on the left of a ~ operator 
#' and the explanatory variables on the right, separated by + operators. The argument corresponds
#' to the argument \code{fixed} in function \code{\link[nlme]{lme}}.
#' @param alpha a one-sides linear formula describing the independent variables used to predict 
#' variance in the "alpha" model 
#' @param pop_data a data frame that needs to comprise the variables
#' named on the right of the ~ operator in \code{fixed}, i.e. the explanatory
#' variables, and \code{pop_domains}.
#' @param pop_domains a character string containing the name of a variable that
#' indicates domains in the population data. The variable can be numeric or
#' a factor but needs to be of the same class as the variable named in
#' \code{smp_domains}.
#' @param pop_subdomains a character string containing the name of a variable that 
#' indicates sub-domains in the population data. When this option is specified, a 
#' two-fold nested error model is used. Defaults to \code{NULL}
#' @param smp_data a data frame that needs to comprise all variables named in
#' \code{fixed} and \code{smp_domains}.
#' @param smp_domains a character string containing the name of a variable
#' that indicates domains in the sample data. The variable can be numeric or a
#' factor but needs to be of the same class as the variable named in
#' \code{pop_domains}.
#' #' @param smp_subdomains a character string containing the name of a variable that 
#' indicates sub-domains in the sample data. When this option is specified, a 
#' two-fold nested error model is used. Defaults to \code{NULL}
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
#' must be at least 1. Defaults to 100. 
#' @param seed an integer to set the seed for the random number generator. For
#' the usage of random number generation, see Details. If seed is set to
#' \code{NULL}, seed is chosen randomly. Defaults to \code{123}.
#' @param custom_indicator a list of functions containing the indicators to be
#' calculated additionally. Such functions must depend on the target variable
#' \code{y}, and optional can depend on \code{pop_weights} and the
#' \code{threshold}. Defaults to \code{NULL}.
#' @param na.rm if \code{TRUE}, observations with \code{NA} values are deleted
#' from the population and sample data. For the EBP procedure complete
#' observations are required. Defaults to \code{FALSE}.
#' @param weights a character string containing the name of a variable that
#' indicates weights in the sample data. If a character string is provided
#' a weighted version of the ebp will be used. The variable has to be numeric.
#' Defaults to \code{NULL}.
#' @param pop_weights a character string containing the name of a variable that
#' indicates population weights in the population data. If a character string
#' is provided weighted indicators are estimated using population weights.
#' The variable has to be numeric. Defaults to \code{NULL}.
#' @param aggregate_to a character string containing the name of a variable from
#' population data that indicates the target domain level for which the
#' results are to be displayed. The variable can be numeric or a factor.
#' Defaults to \code{NULL}.
#' @param weights_type a character string. Two different methods for survey
#' weights are available (i) EBP under informative sampling from
#' \cite{Guadarrama et al. (2018)} ("Guadarrama"); (ii) considering survey
#' weights by using the weighting options of \code{\link{nlme}} from
#' \cite{Pinheiro and Bates (2023)} ("nlme"); (iii) considering survey
#' weights by using the weighting options of \code{\link{nlme}} and use these
#' weights also to determine the optimal transformation parameter lambda
#' ("nlme_lambda"). Defaults to \code{"Guadarrama"}.
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
#' @param benchmark_level a character indicating the level at which the
#' benchmarking is performed. This name must be represented in the sample and
#' population data as variable name.
#' @param benchmark_weights the name of variable containing benchmark weights.
#' This is only possible for internal benchmarking and enable users to benchmark
#' with weights differing from the survey weights (Default for weighting for
#' internal benchmarking).
#' @param random_method selects the method used by the plm function to estimate the random effects model.
#' It can take on the values "swar", "amemiya", "walhus" (default), or "nerlove". 
#' @param rescale_weights a logical indicating if the sample weights are scaled.
#' If \code{FALSE} (default), the sample weights do not change. When \code{TRUE}
#' , the sample weights are rescaled such that the average weight is 1
#' within each domain.
#' @param Ydump a string specifying the name of a .csv file to save all simulated
#' values of the dependent value, model predictions, and error terms used for
#' point estimation.
#' @param errors a string containing either "normal" or "nonnormal". If normal, error terms 
#' are drawn from a normal distribution. If non-normal, error terms are drawn via a 
#' non-parametric bootstrap. 
#' @param model_parameters a string specifying "fixed" or "variable". If variable is specified,
#' estimates of model parameters beta and sigma will be drawn from their estimated 
#' distribution. Otherwise they are assumed fixed. Defaults to "variable".  
#' @param indicators a list of strings containing outcome indicators that should be calculated. 
#' Defaults to NULL, which selects all indicators.  
#' @return An object of class "ell", "emdi" that provides estimators for
#' regional disaggregated indicators.
#' Several generic functions have methods for the returned object. For a full
#' list and descriptions of the components of objects of class "emdi",
#' see \code{\link{emdiObject}}.
#' @details For Monte-Carlo approximationsrandom number generation is used. 
#' Thus, a seed is set by the
#' argument \code{seed}. \cr \cr
#' The set of predefined indicators includes the mean, median, four further
#' quantiles (10\%, 25\%, 75\% and 90\%), head count ratio, poverty gap, Gini
#' coefficient and the quintile share ratio. \cr \cr
#' Since the sample observations often cannot be identified in practical
#' applications, a modified approach by Guadarrama et al. (2016) called census
#' EBP is implemented for the point estimation. For the MSE estimation, the
#' bootstrap sample is not extracted from the superpopulation, but generated by
#' the estimated model parameters. The lower the ratio between the sample and
#' the population size, the closer are the results to the proposed approach by
#' @export

ell <- function(fixed,
                alpha = NULL, 
                pop_data,
                pop_domains,
                pop_subdomains = NULL, 
                smp_data,
                smp_domains,
                smp_subdomains = NULL, 
                L = 100,
                threshold = NULL,
                transformation = "log.shift",
                interval = "default",
                seed = 123,
                custom_indicator = NULL,
                na.rm = FALSE,
                weights = NULL,
                weights_type = NULL, 
                pop_weights = NULL,
                aggregate_to = NULL,
                benchmark = NULL,
                benchmark_type = "ratio",
                benchmark_level = NULL,
                benchmark_weights = NULL,
                rescale_weights = FALSE,
                random_method = "walhus", 
                Ydump = NULL, 
                errors = "normal",
                model_parameters = "variable", 
                indicators = NULL
) {
  
  start.time <- Sys.time()
  ell_check1(
    fixed = fixed, pop_data = pop_data, pop_domains = pop_domains,
    smp_data = smp_data, smp_domains = smp_domains, L = L
  )
  
  
  ell_check2(
    threshold = threshold, transformation = transformation,
    interval = interval, MSE = F, boot_type = "parametric", B = 50, L= L, 
    custom_indicator = custom_indicator, cpus = 1, seed = seed,
    na.rm = na.rm, weights = weights, pop_weights = pop_weights,
    weights_type = weights_type, benchmark = benchmark,
    benchmark_type = benchmark_type, benchmark_level = benchmark_level,
    benchmark_weights = benchmark_weights,MSE_pop_weights=NULL 
  )
  
  if (!is.null(alpha) & transformation!="log" & transformation!="log.shift" & transformation!="no") {
    stop(strwrap(prefix = " ", initial = "",
                 "Must use log or log shift transformation when using alpha model."))
  }
  
  # Save function call ---------------------------------------------------------
  
  call <- match.call()
  if (inherits(call$fixed, "name")) {
    call$fixed <- fixed
  }
  # Data manipulation and notational framework ---------------------------------
  set.seed(seed)

  if (is.null(benchmark_weights) & !is.null(weights)) {
    benchmark_weights <- weights
  }
  
  # The function framework_ell can be found in script framework_ell.R
  framework <- framework_ell(
    pop_data = pop_data,
    pop_domains = pop_domains,
    pop_subdomains = pop_subdomains, 
    smp_data = smp_data,
    smp_domains = smp_domains,
    smp_subdomains = smp_subdomains, 
    aggregate_to = aggregate_to,
    custom_indicator = custom_indicator,
    fixed = fixed,
    alpha = alpha, 
    threshold = threshold,
    na.rm = na.rm,
    weights = weights,
    pop_weights = pop_weights,
    benchmark_level = benchmark_level,
    benchmark_weights = benchmark_weights,
    rescale_weights = rescale_weights,
    errors = errors,
    indicators = indicators,
    random_method = random_method, 
    model_parameters = model_parameters
  )
  
  
  
  # Point Estimation -----------------------------------------------------------
  # The function point_estim can be found in script point_estimation_ell.R
  point_estim <- point_estim_ell(
    framework = framework,
    fixed = fixed,
    alpha = alpha, 
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
  
  
  

    ell_out <- list(
      ind = point_estim$ind,
      var = point_estim$var,
      transform_param = point_estim[c(
        "optimal_lambda",
        "shift_par"
      )],
      model = point_estim$model,
      alpha_model = point_estim$alpha_model, 
      model_par = point_estim$model_par,
      framework = framework[c(
        "N_dom_unobs",
        "N_dom_smp",
        "N_smp",
        "N_pop",
        "smp_domains",
        "smp_data",
        "smp_domains_vec",
        "pop_domains_vec",
        "n_smp",
        "n_pop",
        "response"
      )],
      transformation = transformation,
      method = "reml",
      fixed = fixed,
      alpha = alpha, 
      call = call,
      successful_bootstraps = NULL
    )
  
  
  
  end.time <- Sys.time()
  print(round(end.time - start.time,2))
  
  
  class(ell_out) <- c("ell","povmap")
  return(ell_out)
}
