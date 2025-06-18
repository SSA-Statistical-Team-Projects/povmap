# ========================================================
# Non-Normal Empirical Bayes Implementation (Version 6)
# ========================================================
# This implementation extends the EBP approach by replacing the single normal
# distribution assumption for area-level random effects with a normal mixture model,
# allowing for better handling of non-normal distributions.
# 
#' Main wrapper function for NNEBP
#'
#' @param fixed formula for the fixed effects model
#' @param pop_data population data frame
#' @param pop_domains name of domain variable in population data
#' @param smp_data sample data frame
#' @param smp_domains name of domain variable in sample data
#' @param L number of Monte Carlo simulations
#' @param threshold poverty threshold (numeric or function)
#' @param transformation transformation type
#' @param interval interval for transformation parameter optimization
#' @param MSE logical, whether to compute MSE
#' @param B number of bootstrap replicates for MSE estimation
#' @param mse_method method for MSE estimation: "variance" or "superpopulation"
#' @param seed random seed
#' @param boot_type bootstrap type
#' @param na.rm logical, whether to remove NA values
#' @param custom_indicator list of custom indicator functions
#' @param weights name of weight variable in sample data
#' @param parallel_mode parallelization mode
#' @param cpus number of CPUs for parallel processing
#' @param Ydump optional path to save household-level simulated values
#' @param Pdump optional path to save area-level simulated values and indicators
#' @param ... additional arguments
#'
#' @return object of class "nnebp" and "emdi"
nnebp <- function(fixed, pop_data, pop_domains, smp_data, smp_domains,
                  L = 50, threshold = NULL, transformation = "box.cox",
                  interval = "default", MSE = FALSE, B = 50, seed = 123,
                  mse_method = "superpopulation", boot_type = "parametric",
                  na.rm = FALSE, custom_indicator = NULL, weights = NULL,
                  parallel_mode = "none", cpus = 1, Ydump = NULL, Pdump = NULL, ...) {
  
  # Record the call for later reference
  call <- match.call()
  
  # Set random seed
  if (!is.null(seed)) set.seed(seed)
  
  # Get response variable name
  response_var <- all.vars(fixed)[1]
  
  # Handle NA values if requested
  if (na.rm) {
    na_idx_smp <- stats::complete.cases(smp_data)
    smp_data <- smp_data[na_idx_smp, ]
    
    # Also handle NAs in population data if needed
    pop_vars <- all.vars(fixed)[-1]  # Exclude response variable
    pop_vars <- pop_vars[pop_vars %in% names(pop_data)]
    na_idx_pop <- stats::complete.cases(pop_data[, c(pop_domains, pop_vars)])
    pop_data <- pop_data[na_idx_pop, ]
  }
  
  # Default threshold - ensure it matches standard EBP approach
  if (is.null(threshold)) {
    cat("Setting default threshold to 60% of median of", response_var, "\n")
    threshold <- function(y) 0.6 * stats::median(y)
    
    # Calculate the actual threshold value for debugging
    threshold_value <- 0.6 * stats::median(smp_data[[response_var]])
    cat("Calculated threshold value:", threshold_value, "\n")
  }
  
  # Create framework that matches the EBP approach
  framework <- list(
    smp_data = smp_data,
    smp_domains = smp_domains,
    pop_data = pop_data,
    pop_domains = pop_domains,
    pop_domains_vec = as.character(pop_data[[pop_domains]]),
    threshold = threshold,
    N_smp = nrow(smp_data),
    N_pop = nrow(pop_data),
    weights = weights,
    model_parameters = "variable",
    custom_indicator = custom_indicator,
    random_method = "amemiya",  # Use a stable method
    fixed = fixed  # Store the formula for later use
  )
  
  # Point estimation
  cat("Running point estimation with NNEBP...\n")
  point_est <- point_estim_nnebp(
    framework = framework,
    fixed = fixed,
    transformation = transformation,
    interval = interval,
    L = L,
    Ydump = Ydump,
    Pdump = Pdump
  )
  
  # MSE estimation if requested
  mse_est <- NULL
  if (MSE) {
    cat("Running MSE estimation with parametric bootstrap...\n")
    cat("Using MSE method:", mse_method, "\n")
    mse_result <- parametric_bootstrap_nnebp(
      framework = framework,
      fixed = fixed,
      transformation = transformation,
      B = B,
      L = L,
      mse_method = mse_method
    )
    mse_est <- mse_result$mse_estimates
  }
  
  # Create result object that matches standard EBP format
  result <- list(
    call = call,
    fixed = fixed,
    framework = framework,
    ind = point_est$ind,
    MSE = mse_est,
    method = "reml",
    model = list(
      model_par = point_est$model_par,
      fit = point_est$model
    ),
    transformation = list(
      type = transformation,
      optimal_lambda = point_est$optimal_lambda,
      backtransformation = TRUE
    ),
    transform_param = list(
      optimal_lambda = point_est$optimal_lambda,
      shift_par = point_est$shift_par
    )
  )
  
  # Set class and return
  class(result) <- c("nnebp", "emdi")
  return(result)
}

#' Point estimation for Non-Normal EB
#'
#' @param framework framework list with data and configuration
#' @param fixed formula for fixed effects
#' @param transformation transformation type
#' @param interval interval for transformation parameter
#' @param L number of Monte Carlo simulations
#' @param Ydump optional path to save household-level simulated values
#' @param Pdump optional path to save area-level simulated values and indicators
#'
#' @return list with point estimates and model information
point_estim_nnebp <- function(framework, fixed, transformation = "box.cox", 
                              interval = "default", L = 50, Ydump = NULL, Pdump = NULL) {
  cat("Starting point estimation...\n")
  
  # Get response variable name
  response_var <- all.vars(fixed)[1]
  cat("Response variable:", response_var, "\n")
  
  # Ensure alignment with standard EBP by using the same optimal parameter approach
  cat("Obtaining optimal transformation parameter...\n")
  optimal_lambda <- optimal_parameter(
    generic_opt = generic_opt,  # assumed to be defined in povmap
    fixed = fixed,
    smp_data = framework$smp_data,
    smp_domains = framework$smp_domains,
    transformation = transformation,
    interval = interval,
    framework = framework
  )
  cat("Optimal lambda:", optimal_lambda, "\n")
  
  # Transform the sample data using the optimal lambda
  cat("Transforming data...\n")
  transformation_par <- data_transformation(
    fixed = fixed,
    smp_data = framework$smp_data,
    transformation = transformation,
    lambda = optimal_lambda
  )
  shift_par <- transformation_par$shift
  cat("Shift parameter:", shift_par, "\n")
  
  # Prepare arguments for fitting the random effects model via plm
  args <- list(
    formula = fixed,
    data = transformation_par$transformed_data,
    model = "random",
    index = framework$smp_domains
  )
  
  # Force the use of a known stable random effects method ("amemiya")
  args$random.method <- "amemiya"
  
  # Fit the random effects model; stop if convergence fails
  cat("Fitting random effects model...\n")
  re_model <- tryCatch({
    do.call(plm::plm, args)
  }, error = function(e) {
    stop("Model did not converge using plm::plm. Error message: ", e$message)
  })
  
  # Obtain model parameters including the estimated normal mixture components
  cat("Extracting model parameters...\n")
  model_par <- model_par_nnebp(
    framework = framework,
    re_model = re_model,
    fixed = fixed,
    transformation_par = transformation_par
  )
  
  # Print normal mixture parameters for debugging
  cat("Normal mixture model parameters:\n")
  cat("Mixture components:", length(model_par$normal_mixture$lambda), "\n")
  cat("Mixing proportions:", model_par$normal_mixture$lambda, "\n")
  cat("Component means:", model_par$normal_mixture$mu, "\n")
  cat("Component std devs:", model_par$normal_mixture$sigma, "\n")
  
  # Run Monte Carlo simulation to generate pseudo responses
  cat("Running Monte Carlo simulation with", L, "iterations...\n")
  indicator_prediction <- monte_carlo_nnebp(
    framework = framework,
    L = L,
    model_par = model_par,
    seed = 123,
    Ydump = Ydump,
    Pdump = Pdump,
    transformation = transformation,
    lambda = optimal_lambda,
    shift = shift_par
  )
  
  # Apply back-transformation to simulated responses
  cat("Back-transforming simulated values...\n")
  simulated_y <- lapply(indicator_prediction$simulated_values, function(y) {
    back_transformation(
      y = y,
      transformation = transformation,
      lambda = optimal_lambda,
      shift = shift_par,
      framework = framework,
      fixed = fixed
    )
  })
  
  # Calculate indicators for each simulation replication
  cat("Calculating indicators...\n")
  if (is.function(framework$threshold)) {
    # If threshold is a function, print what we'll be using
    sample_y <- simulated_y[[1]]
    threshold_value <- framework$threshold(sample_y)
    cat("Using threshold function, example value:", threshold_value, "\n")
  } else {
    cat("Using fixed threshold value:", framework$threshold, "\n")
  }
  
  simulated_indicators <- lapply(simulated_y, function(y) {
    calculate_indicators(
      y = y,
      framework = framework,
      threshold = framework$threshold
    )
  })
  
  # Average the indicators elementwise across replications
  indicator_names <- names(simulated_indicators[[1]])
  point_estimates <- lapply(indicator_names, function(nm) {
    replicate_matrix <- do.call(rbind, lapply(simulated_indicators, function(ind) ind[[nm]]))
    colMeans(replicate_matrix, na.rm = TRUE)
  })
  names(point_estimates) <- indicator_names
  
  # Calculate the variance of indicators across replications
  var_estimates <- lapply(indicator_names, function(nm) {
    replicate_matrix <- do.call(rbind, lapply(simulated_indicators, function(ind) ind[[nm]]))
    apply(replicate_matrix, 2, var, na.rm = TRUE)
  })
  names(var_estimates) <- indicator_names
  
  # Debug output
  cat("Summary of indicator estimates:\n")
  for (ind in indicator_names) {
    cat(ind, "- Mean:", mean(point_estimates[[ind]], na.rm = TRUE), 
        "Min:", min(point_estimates[[ind]], na.rm = TRUE), 
        "Max:", max(point_estimates[[ind]], na.rm = TRUE), "\n")
  }
  
  list(
    ind = point_estimates,
    var = var_estimates,
    optimal_lambda = optimal_lambda,
    shift_par = shift_par,
    model = re_model,
    model_par = model_par
  )
}

#' Monte Carlo Simulation with Conditional Normal Mixture
#'
#' @param framework framework list with data and configuration
#' @param L number of Monte Carlo simulations
#' @param model_par model parameters
#' @param seed random seed
#' @param Ydump optional path to save household-level simulated values
#' @param Pdump optional path to save area-level simulated values and indicators
#' @param transformation transformation type
#' @param lambda transformation parameter
#' @param shift shift parameter
#'
#' @return list with simulated values
monte_carlo_nnebp <- function(framework, L, model_par, seed = 42, Ydump = NULL, Pdump = NULL,
                              transformation = NULL, lambda = NULL, shift = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  # Initialize Ydump file if specified
  if (!is.null(Ydump)) {
    Ydumpdf <- data.frame(matrix(ncol = 6, nrow = 0))
    colnames(Ydumpdf) <- c("L", "Domain", "Simulated_Y", "XBetahat", "eta", "epsilon")
    write.csv(Ydumpdf, Ydump, row.names = FALSE)
  }
  
  # Initialize Pdump file if specified
  if (!is.null(Pdump)) {
    # Define columns for area-level data
    Pdumpdf <- data.frame(matrix(ncol = 13, nrow = 0))
    colnames(Pdumpdf) <- c("L", "Domain", "Area_Effect", "Mean", "Head_Count", 
                           "Poverty_Gap", "Gini", "Median", "Quantile_10", 
                           "Quantile_25", "Quantile_75", "Quantile_90", "Quintile_Share")
    write.csv(Pdumpdf, Pdump, row.names = FALSE)
  }
  
  # Create design matrix for the population data
  covariates_terms <- attr(terms(framework$fixed), "term.labels")
  new_formula <- as.formula(paste("~", paste(covariates_terms, collapse = " + ")))
  X_pop <- model.matrix(new_formula, framework$pop_data)
  mu_pop <- as.vector(X_pop %*% model_par$betas)
  
  # Get unique areas for area-level output
  unique_areas <- unique(as.character(framework$pop_domains_vec))
  
  # Initialize simulations list
  simulations <- list()
  
  # Get population weights if available
  if (!is.null(framework$pop_weights)) {
    pop_weights_vec <- framework$pop_data[[framework$pop_weights]]
  } else {
    pop_weights_vec <- rep(1, length(framework$pop_domains_vec))
  }
  
  for (l in 1:L) {
    # Generate errors
    errors <- errors_gen_nnebp(
      framework = framework, 
      model_par = model_par, 
      unit_level = TRUE, 
      seed = NULL
    )
    
    # Extract area effects and idiosyncratic errors
    area_effects_pop <- errors$area_effects[framework$pop_domains_vec]
    idio_errors_pop <- errors$idio_errors
    
    # Combine all components for the population
    y_sim <- mu_pop + area_effects_pop + idio_errors_pop
    
    # Attach domain labels
    attributes(y_sim)$domains <- framework$pop_domains_vec
    
    # Store simulation
    simulations[[l]] <- y_sim
    
    # Output household-level data if Ydump is specified
    if (!is.null(Ydump)) {
      Ydumpdf <- data.frame(rep(l, length(y_sim)), 
                            framework$pop_domains_vec,
                            y_sim,
                            mu_pop,
                            area_effects_pop,
                            idio_errors_pop)
      write.table(Ydumpdf, file = Ydump, row.names = FALSE, 
                  append = TRUE, col.names = FALSE, sep = ",")
    }
    
    # Calculate and output area-level indicators if Pdump is specified
    if (!is.null(Pdump)) {
      # Apply back-transformation if transformation parameters provided
      if (!is.null(transformation) && !is.null(lambda)) {
        y_back <- back_transformation(
          y = y_sim,
          transformation = transformation,
          lambda = lambda,
          shift = shift,
          framework = framework,
          fixed = framework$fixed
        )
      } else {
        y_back <- y_sim  # Use untransformed values if parameters not provided
      }
      
      # Calculate indicators by domain for this simulation
      area_indicators <- calculate_indicators(
        y = y_back,
        framework = framework,
        threshold = framework$threshold
      )
      
      # Extract area effects for each domain
      domain_effects <- sapply(unique_areas, function(domain) {
        mean(errors$area_effects[names(errors$area_effects) == domain])
      })
      
      # Create data for each area
      for (i in seq_along(unique_areas)) {
        domain <- unique_areas[i]
        area_effect <- domain_effects[i]
        
        # Extract indicators for this domain
        mean_val <- area_indicators$Mean[i]
        hc_val <- area_indicators$Head_Count[i]
        pg_val <- area_indicators$Poverty_Gap[i]
        gini_val <- area_indicators$Gini[i]
        median_val <- area_indicators$Median[i]
        q10_val <- area_indicators$Quantile_10[i]
        q25_val <- area_indicators$Quantile_25[i]
        q75_val <- area_indicators$Quantile_75[i]
        q90_val <- area_indicators$Quantile_90[i]
        qs_val <- area_indicators$Quintile_Share[i]
        
        # Create row for this area
        Pdumpdf <- data.frame(
          l, domain, area_effect, mean_val, hc_val, pg_val, gini_val,
          median_val, q10_val, q25_val, q75_val, q90_val, qs_val
        )
        
        # Append to Pdump file
        write.table(Pdumpdf, file = Pdump, row.names = FALSE,
                    append = TRUE, col.names = FALSE, sep = ",")
      }
    }
  }
  
  list(simulated_values = simulations)
}

#' Model parameter estimation for NNEBP
#'
#' @param framework framework list with data and configuration
#' @param re_model fitted random effects model
#' @param fixed formula for fixed effects
#' @param transformation_par transformation parameters
#'
#' @return list of model parameters
model_par_nnebp <- function(framework, re_model, fixed, transformation_par) {
  require(mixtools)
  
  if (is.null(framework$smp_data) || is.null(framework$smp_domains)) {
    stop("Sample data and domains are required in framework.")
  }
  
  # Extract model parameters from the fitted model
  betas <- re_model$coefficients 
  sigmae2est <- re_model$ercomp$sigma2[1]  
  sigmau2est <- re_model$ercomp$sigma2[2]  
  
  # Handle potential NULL model parameters field
  varFix <- if (!is.null(framework$model_parameters) && framework$model_parameters != "fixed") {
    re_model$vcov 
  } else {
    NULL
  }
  
  # Calculate residuals and handle weights (defaulting to 1 if not specified)
  residuals <- as.vector(re_model$model[, 1] - predict(re_model))
  weights <- if (is.null(framework$weights)) {
    rep(1, length(residuals)) 
  } else {
    framework$smp_data[, framework$weights]
  }
  
  # Aggregate weighted residuals by domain
  agg <- aggregate_weighted_mean(
    df = residuals,
    by = list(area = framework$smp_data[, framework$smp_domains]),
    w = weights
  )
  mean_residuals <- setNames(agg[, 2], agg[, 1])
  mean_residuals <- mean_residuals[!is.na(mean_residuals)]
  
  # Estimate a normal mixture model from the aggregated residuals
  # Default to a single component if unable to fit a mixture
  tryCatch({
    nm_result <- estimate_normal_mixture(mean_residuals)
    cat("Fitted normal mixture with", length(nm_result$lambda), "components\n")
  }, error = function(e) {
    cat("Error fitting normal mixture:", e$message, "\n")
    cat("Falling back to single normal component\n")
    nm_result <<- list(
      lambda = 1,
      mu = mean(mean_residuals),
      sigma = sd(mean_residuals),
      convergence = TRUE,
      k = 1
    )
  })
  
  # Compute conditional parameters for each domain
  conditional_params <- compute_conditional_params(
    nm_result = nm_result,
    mean_residuals = mean_residuals,
    sigmae2est = sigmae2est,
    framework = framework
  )
  
  # Compute fixed effect predictions
  X_smp <- model.matrix(re_model)
  mu_fixed <- as.vector(X_smp %*% betas)
  
  list(
    betas = betas,
    sigmae2est = sigmae2est,
    sigmau2est = sigmau2est,
    varFix = varFix,
    residuals = residuals,
    mean_residuals = mean_residuals,
    mu_fixed = mu_fixed,
    normal_mixture = nm_result,
    conditional_params = conditional_params
  )
}

#' Generate errors from normal mixture model
#'
#' @param framework framework list with data and configuration
#' @param model_par model parameters
#' @param areas specific areas to generate errors for (NULL = all)
#' @param unit_level logical, whether to generate unit-level errors
#' @param seed random seed
#'
#' @return list with area effects and idiosyncratic errors
errors_gen_nnebp <- function(framework, model_par, areas = NULL, unit_level = FALSE, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  # Determine which areas to generate errors for
  if (is.null(areas)) {
    # Use all areas
    unique_areas <- unique(as.character(framework$pop_data[, framework$pop_domains]))
  } else {
    # Use specified areas
    unique_areas <- areas
  }
  
  # Generate area effects from the normal mixture
  area_effects <- sapply(unique_areas, function(area) {
    # Get conditional parameters for this area
    params <- model_par$conditional_params[[area]]
    
    # If this is an out-of-sample area, use unconditional distribution
    if (is.null(params)) {
      # For out-of-sample areas, use unconditional mixture
      nm <- model_par$normal_mixture
      
      # Sample a component based on mixing weights
      comp_index <- sample(seq_along(nm$lambda), size = 1, prob = nm$lambda)
      
      # Draw from the selected component
      rnorm(1, mean = nm$mu[comp_index], sd = nm$sigma[comp_index])
    } else {
      # For in-sample areas, use conditional mixture
      # Sample a component based on conditional weights
      comp_index <- sample(seq_along(params$weights), size = 1, prob = params$weights)
      
      # Draw from the selected component
      rnorm(1, mean = params$means[comp_index], sd = sqrt(params$variances[comp_index]))
    }
  })
  names(area_effects) <- unique_areas
  
  # Generate idiosyncratic errors if requested
  idio_errors <- NULL
  if (unit_level) {
    # For unit level, we need one error per population unit
    # Sample from the empirical distribution of residuals
    idio_errors <- sample(
      model_par$residuals, 
      size = nrow(framework$pop_data), 
      replace = TRUE
    )
  }
  
  list(
    area_effects = area_effects,
    idio_errors = idio_errors
  )
}

#' Helper function to generate superpopulation
#'
#' @param framework framework list with data and configuration
#' @param model_par model parameters
#' @param seed random seed
#'
#' @return superpopulation data
generate_superpopulation <- function(framework, model_par, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  cat("Generating superpopulation...\n")
  
  # Extract population dimensions
  n_pop <- nrow(framework$pop_data)
  
  # Create design matrix for population
  covariates_terms <- attr(terms(framework$fixed), "term.labels")
  new_formula <- as.formula(paste("~", paste(covariates_terms, collapse = " + ")))
  X_pop <- model.matrix(new_formula, framework$pop_data)
  
  # Calculate fixed effects
  mu_pop <- as.vector(X_pop %*% model_par$betas)
  
  # Generate area effects for each unit based on the normal mixture model
  pop_domains <- as.character(framework$pop_data[, framework$pop_domains])
  unique_domains <- unique(pop_domains)
  
  # Generate one area effect per domain from the unconditional mixture
  nm <- model_par$normal_mixture
  area_effects <- sapply(unique_domains, function(d) {
    # Sample component based on mixing weights
    comp_idx <- sample(1:length(nm$lambda), 1, prob = nm$lambda)
    # Draw from selected component
    rnorm(1, mean = nm$mu[comp_idx], sd = nm$sigma[comp_idx])
  })
  names(area_effects) <- unique_domains
  
  # Expand area effects to all population units
  area_effects_pop <- area_effects[pop_domains]
  
  # Generate idiosyncratic errors for all population units
  idio_errors_pop <- rnorm(n_pop, mean = 0, sd = sqrt(model_par$sigmae2est))
  
  # Combine all components
  y_sim <- mu_pop + area_effects_pop + idio_errors_pop
  
  # Return the superpopulation
  superpop <- framework$pop_data
  superpop$y_sim <- y_sim
  
  cat("Superpopulation generated with", n_pop, "units\n")
  return(list(
    superpop = superpop,
    area_effects = area_effects,
    mu_pop = mu_pop
  ))
}

#' Parametric Bootstrap for MSE Estimation with NNEBP
#'
#' @param framework framework list with data and configuration
#' @param fixed formula for fixed effects
#' @param transformation transformation type
#' @param B number of bootstrap replicates
#' @param L number of Monte Carlo simulations
#' @param mse_method method for MSE estimation ("variance" or "superpopulation")
#'
#' @return list with MSE estimates
parametric_bootstrap_nnebp <- function(framework, fixed, transformation = "box.cox", 
                                       B = 100, L = 50, mse_method = "superpopulation") {
  cat("Starting parametric bootstrap for MSE estimation...\n")
  cat("Using MSE method:", mse_method, "\n")
  
  # Step 1: Obtain transformation parameters and fit the original model
  cat("Obtaining transformation parameters...\n")
  optimal_lambda <- optimal_parameter(
    generic_opt = generic_opt,
    fixed = fixed,
    smp_data = framework$smp_data,
    smp_domains = framework$smp_domains,
    transformation = transformation,
    interval = "default",
    framework = framework
  )
  
  transformation_par <- data_transformation(
    fixed = fixed,
    smp_data = framework$smp_data,
    transformation = transformation,
    lambda = optimal_lambda
  )
  shift_par <- transformation_par$shift
  
  # Build arguments for plm
  args <- list(
    formula = fixed,
    data = transformation_par$transformed_data,
    model = "random",
    index = framework$smp_domains
  )
  valid_methods <- c("swar", "walhus", "amemiya", "nerlove")
  if (!is.null(framework$random_method) && framework$random_method %in% valid_methods) {
    args$random.method <- framework$random_method
  }
  
  cat("Fitting original model...\n")
  re_model <- do.call(plm::plm, args)
  
  model_par <- model_par_nnebp(
    framework = framework,
    re_model = re_model,
    fixed = fixed,
    transformation_par = transformation_par
  )
  
  # Compute original point estimates
  cat("Computing original point estimates...\n")
  original_point <- point_estim_nnebp(framework, fixed, transformation, interval = "default", L = L)
  original_indicators <- original_point$ind
  
  # Branch based on MSE method
  if (mse_method == "superpopulation") {
    # Superpopulation approach for MSE
    cat("Using superpopulation approach for MSE...\n")
    
    # Step 2: Generate superpopulation - "true" population based on the fitted model
    cat("Generating superpopulation...\n")
    superpop_result <- generate_superpopulation(framework, model_par, seed = 123)
    superpop <- superpop_result$superpop
    
    # Extract response variable name
    response_var <- all.vars(fixed)[1]
    
    # Create a framework for the superpopulation
    superpop_framework <- framework
    superpop_framework$pop_data <- superpop
    
    # Step 3: Calculate "true" indicator values from the superpopulation
    cat("Calculating true indicator values from superpopulation...\n")
    
    # For true indicator values, we use simulated values
    y_true <- superpop$y_sim
    
    # Apply backtransformation if needed
    if (transformation == "box.cox") {
      if (abs(optimal_lambda) < 1e-6) {
        y_true <- exp(y_true)
      } else {
        y_true <- ((optimal_lambda * y_true) + 1)^(1/optimal_lambda)
      }
    } else if (transformation == "log") {
      y_true <- exp(y_true)
    }
    
    # Set domain attribute for indicator calculation
    attr(y_true, "domains") <- as.character(superpop[, framework$pop_domains])
    
    # Calculate true indicators
    true_indicators <- calculate_indicators(
      y = y_true,
      framework = superpop_framework,
      threshold = framework$threshold
    )
    
    # Initialize arrays to store squared errors
    indicator_names <- names(true_indicators)
    n_domains <- length(true_indicators[[1]])
    squared_errors <- list()
    for (ind_name in indicator_names) {
      squared_errors[[ind_name]] <- matrix(0, nrow = B, ncol = n_domains)
    }
    
    # Step 4: For each bootstrap replicate, generate a sample from the superpopulation
    cat("Running", B, "bootstrap replications...\n")
    for (b in seq_len(B)) {
      cat("Bootstrap replication", b, "of", B, "\n")
      
      # Create a bootstrap sample from the superpopulation
      # This mimics how the original sample was drawn from the population
      unique_domains <- unique(superpop[, framework$pop_domains])
      bootstrap_sample <- data.frame()
      
      # For each domain, sample with the same sampling rate as in the original sample
      for (domain in unique_domains) {
        # Get number of units in original sample for this domain
        n_orig <- sum(framework$smp_data[, framework$smp_domains] == domain)
        
        # Get units in superpopulation for this domain
        domain_units <- superpop[superpop[, framework$pop_domains] == domain, ]
        
        # Sample n_orig units with replacement
        if (nrow(domain_units) > 0 && n_orig > 0) {
          domain_sample <- domain_units[sample(1:nrow(domain_units), n_orig, replace = TRUE), ]
          bootstrap_sample <- rbind(bootstrap_sample, domain_sample)
        }
      }
      
      # Skip this iteration if bootstrap sample is empty
      if (nrow(bootstrap_sample) == 0) {
        cat("Warning: Empty bootstrap sample in replicate", b, "- skipping\n")
        next
      }
      
      # Replace original y-variable with simulated y
      bootstrap_sample[[response_var]] <- bootstrap_sample$y_sim
      
      # Update the framework for the bootstrap replicate
      framework_boot <- framework
      framework_boot$smp_data <- bootstrap_sample
      
      # Transform the bootstrap sample
      transformation_par_boot <- tryCatch({
        data_transformation(
          fixed = fixed,
          smp_data = bootstrap_sample,
          transformation = transformation,
          lambda = optimal_lambda
        )
      }, error = function(e) {
        cat("Warning: Data transformation failed in bootstrap replicate", b, "-", e$message, "\n")
        NULL
      })
      
      # Skip this iteration if transformation failed
      if (is.null(transformation_par_boot)) {
        next
      }
      
      # Fit model to bootstrap sample
      args_boot <- list(
        formula = fixed,
        data = transformation_par_boot$transformed_data,
        model = "random",
        index = framework_boot$smp_domains
      )
      if (!is.null(framework_boot$random_method) && framework_boot$random_method %in% valid_methods) {
        args_boot$random.method <- framework_boot$random_method
      }
      
      # Handle potential convergence issues
      re_model_boot <- tryCatch({
        do.call(plm::plm, args_boot)
      }, error = function(e) {
        cat("Warning: Model fitting failed in bootstrap replicate", b, "-", e$message, "\n")
        NULL
      })
      
      # Skip this iteration if model fitting failed
      if (is.null(re_model_boot)) {
        next
      }
      
      # Extract model parameters from bootstrap model
      model_par_boot <- tryCatch({
        model_par_nnebp(
          framework = framework_boot,
          re_model = re_model_boot,
          fixed = fixed,
          transformation_par = transformation_par_boot
        )
      }, error = function(e) {
        cat("Warning: Parameter extraction failed in bootstrap replicate", b, "-", e$message, "\n")
        NULL
      })
      
      # Skip this iteration if parameter extraction failed
      if (is.null(model_par_boot)) {
        next
      }
      
      # Calculate estimates using the bootstrap model
      bootstrap_point <- tryCatch({
        point_estim_nnebp(framework_boot, fixed, transformation, interval = "default", L = L)
      }, error = function(e) {
        cat("Warning: Point estimation failed in bootstrap replicate", b, "-", e$message, "\n")
        NULL
      })
      
      # If bootstrap point estimation succeeded, calculate squared errors
      if (!is.null(bootstrap_point)) {
        for (ind_name in indicator_names) {
          if (ind_name %in% names(bootstrap_point$ind)) {
            # Calculate squared errors between bootstrap estimates and true values
            squared_errors[[ind_name]][b, ] <- (bootstrap_point$ind[[ind_name]] - true_indicators[[ind_name]])^2
          }
        }
      }
    }
    
    # Step 5: Compute MSE as the mean of squared errors for each indicator and domain
    cat("Computing MSE estimates...\n")
    mse_estimates <- list()
    
    for (ind_name in indicator_names) {
      # Calculate mean of squared errors (MSE)
      mse_estimates[[ind_name]] <- colMeans(squared_errors[[ind_name]], na.rm = TRUE)
      
      # Print summary statistics for MSE
      cat("MSE for", ind_name, "- Mean:", mean(mse_estimates[[ind_name]], na.rm = TRUE), 
          "Min:", min(mse_estimates[[ind_name]], na.rm = TRUE), 
          "Max:", max(mse_estimates[[ind_name]], na.rm = TRUE), "\n")
    }
    
    result <- list(
      original_indicators = original_indicators,
      true_indicators = true_indicators,
      mse_estimates = mse_estimates,
      optimal_lambda = optimal_lambda,
      shift_par = shift_par,
      model = re_model,
      model_par = model_par
    )
    
  } else {
    # Variance-based approach for MSE (original method)
    cat("Using variance-based approach for MSE...\n")
    
    # Initialize list to hold bootstrap replicate indicators
    bootstrap_indicators <- vector("list", B)
    
    # For each bootstrap replicate, generate a pseudo sample and re-run estimation
    cat("Running", B, "bootstrap replications...\n")
    for (b in seq_len(B)) {
      cat("Bootstrap replication", b, "of", B, "\n")
      
      # Generate a pseudo sample from the fitted model
      pseudo_y <- simulate_pseudo_sample(framework, model_par, seed = b)
      
      # Replace the response in the original sample data with the pseudo responses
      pseudo_data <- framework$smp_data
      pseudo_data[, 1] <- pseudo_y
      
      # Update the framework for the bootstrap replicate
      framework_boot <- framework
      framework_boot$smp_data <- pseudo_data
      
      # Re-run the transformation and model fitting
      transformation_par_boot <- tryCatch({
        data_transformation(
          fixed = fixed,
          smp_data = pseudo_data,
          transformation = transformation,
          lambda = optimal_lambda
        )
      }, error = function(e) {
        cat("Warning: Data transformation failed in bootstrap replicate", b, "-", e$message, "\n")
        NULL
      })
      
      # Skip this iteration if transformation failed
      if (is.null(transformation_par_boot)) {
        next
      }
      
      # Fit model to bootstrap sample
      args_boot <- list(
        formula = fixed,
        data = transformation_par_boot$transformed_data,
        model = "random",
        index = framework_boot$smp_domains
      )
      if (!is.null(framework_boot$random_method) && framework_boot$random_method %in% valid_methods) {
        args_boot$random.method <- framework_boot$random_method
      }
      
      # Handle potential convergence issues
      re_model_boot <- tryCatch({
        do.call(plm::plm, args_boot)
      }, error = function(e) {
        cat("Warning: Model fitting failed in bootstrap replicate", b, "-", e$message, "\n")
        NULL
      })
      
      # Skip this iteration if model fitting failed
      if (is.null(re_model_boot)) {
        next
      }
      
      # Extract model parameters from bootstrap model
      model_par_boot <- tryCatch({
        model_par_nnebp(
          framework = framework_boot,
          re_model = re_model_boot,
          fixed = fixed,
          transformation_par = transformation_par_boot
        )
      }, error = function(e) {
        cat("Warning: Parameter extraction failed in bootstrap replicate", b, "-", e$message, "\n")
        NULL
      })
      
      # Skip this iteration if parameter extraction failed
      if (is.null(model_par_boot)) {
        next
      }
      
      # Calculate estimates using the bootstrap model
      bootstrap_point <- tryCatch({
        point_estim_nnebp(framework_boot, fixed, transformation, interval = "default", L = L)
      }, error = function(e) {
        cat("Warning: Point estimation failed in bootstrap replicate", b, "-", e$message, "\n")
        NULL
      })
      
      # Store bootstrap estimates if successful
      if (!is.null(bootstrap_point)) {
        bootstrap_indicators[[b]] <- bootstrap_point$ind
      }
    }
    
    # Remove NULL entries from bootstrap_indicators
    bootstrap_indicators <- bootstrap_indicators[!sapply(bootstrap_indicators, is.null)]
    
    # Compute bootstrap MSE estimates for each indicator by domain
    cat("Computing MSE estimates...\n")
    indicator_names <- names(original_indicators)
    mse_estimates <- list()
    
    for (ind_name in indicator_names) {
      # Create matrix of bootstrap estimates
      bootstrap_matrix <- do.call(rbind, lapply(bootstrap_indicators, function(ind_list) {
        ind_list[[ind_name]]
      }))
      
      # Calculate variance across bootstrap replicates
      mse_estimates[[ind_name]] <- apply(bootstrap_matrix, 2, var, na.rm = TRUE)
      
      # Print summary statistics for MSE
      cat("MSE for", ind_name, "- Mean:", mean(mse_estimates[[ind_name]], na.rm = TRUE), 
          "Min:", min(mse_estimates[[ind_name]], na.rm = TRUE), 
          "Max:", max(mse_estimates[[ind_name]], na.rm = TRUE), "\n")
    }
    
    result <- list(
      original_indicators = original_indicators,
      bootstrap_indicators = bootstrap_indicators,
      mse_estimates = mse_estimates,
      optimal_lambda = optimal_lambda,
      shift_par = shift_par,
      model = re_model,
      model_par = model_par
    )
  }
  
  return(result)
}

#' Generate pseudo sample for bootstrap
#'
#' @param framework framework list with data and configuration
#' @param model_par model parameters
#' @param seed random seed
#'
#' @return simulated response values
simulate_pseudo_sample <- function(framework, model_par, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  # Generate sample-level errors 
  sample_areas <- as.character(framework$smp_data[, framework$smp_domains])
  
  # Use errors_gen_nnebp to generate the errors
  errors <- errors_gen_nnebp(
    framework = framework,
    model_par = model_par,
    areas = unique(sample_areas),
    unit_level = FALSE,  # We'll handle idiosyncratic errors separately
    seed = NULL  # Already set at function level
  )
  
  # Extract area effects for sample observations
  area_effects_sim <- errors$area_effects[sample_areas]
  
  # Draw idiosyncratic errors for sample
  idio_errors <- sample(model_par$residuals, size = nrow(framework$smp_data), replace = TRUE)
  
  # Combine fixed effects, area effects, and idiosyncratic errors
  y_sim <- model_par$mu_fixed + area_effects_sim + idio_errors
  
  # Attach domain attribute for consistency
  attr(y_sim, "domains") <- framework$pop_domains_vec[1:length(y_sim)]
  
  # Return the result
  cat("Generated pseudo sample with mean:", mean(y_sim, na.rm = TRUE), "\n")
  y_sim
}

#' Estimate Normal Mixture Distribution
#'
#' @param area_effects vector of area-level effects
#' @param max_components maximum number of mixture components
#' @param min_variance minimum variance for components
#'
#' @return normal mixture model parameters
estimate_normal_mixture <- function(area_effects, max_components = 2, min_variance = 0.001) {
  require(mixtools)
  
  if (!is.numeric(area_effects)) stop("area_effects must be numeric")
  if (max_components < 1) stop("max_components must be positive")
  if (min_variance <= 0) stop("min_variance must be positive")
  
  area_effects <- area_effects[!is.na(area_effects)]
  if (length(area_effects) < 2) {
    return(list(
      lambda = 1,
      mu = mean(area_effects, na.rm = TRUE),
      sigma = sd(area_effects, na.rm = TRUE),
      convergence = TRUE,
      k = 1
    ))
  }
  
  get_init_values <- function(k) {
    if (k == 1) {
      return(list(
        mu = mean(area_effects),
        sigma = max(sd(area_effects), sqrt(min_variance)),
        lambda = 1
      ))
    }
    
    data_range <- range(area_effects)
    spread <- diff(data_range)
    mu_init <- seq(data_range[1] + 0.1 * spread, data_range[2] - 0.1 * spread, length.out = k)
    sigma_init <- rep(max(sd(area_effects) / sqrt(k), sqrt(min_variance)), k)
    
    list(
      mu = mu_init,
      sigma = sigma_init,
      lambda = rep(1/k, k)
    )
  }
  
  mixture_fits <- list()
  n <- length(area_effects)
  
  # Fit k = 1 model
  mixture_fits[[1]] <- list(
    lambda = 1,
    mu = mean(area_effects),
    sigma = max(sd(area_effects), sqrt(min_variance)),
    convergence = TRUE,
    k = 1,
    bic = -2 * sum(dnorm(area_effects, mean(area_effects), sd(area_effects), log = TRUE)) + log(n)
  )
  
  # Fit models for k > 1
  for (k in 2:max_components) {
    init <- get_init_values(k)
    
    fit <- tryCatch({
      em_fit <- normalmixEM(
        area_effects,
        lambda = init$lambda,
        mu = init$mu,
        sigma = init$sigma,
        k = k,
        maxit = 1000,
        epsilon = 1e-6
      )
      
      if (any(em_fit$sigma^2 < min_variance) || any(is.na(em_fit$lambda)) || any(em_fit$lambda < 0.01)) {
        return(NULL)
      }
      
      comp_dens <- sapply(1:k, function(j) {
        em_fit$lambda[j] * dnorm(area_effects, em_fit$mu[j], em_fit$sigma[j])
      })
      ll <- sum(log(rowSums(comp_dens)))
      bic <- -2 * ll + (3 * k - 1) * log(n)
      
      em_fit$convergence <- TRUE
      em_fit$k <- k
      em_fit$bic <- bic
      em_fit
    }, warning = function(w) NULL, error = function(e) NULL)
    
    if (!is.null(fit)) mixture_fits[[k]] <- fit
  }
  
  bic_scores <- sapply(mixture_fits, function(x) { if (is.null(x)) return(Inf) else x$bic })
  best_model <- mixture_fits[[which.min(bic_scores)]]
  best_model$fits_attempted <- mixture_fits
  best_model$bic_scores <- bic_scores
  
  best_model
}

#' Compute Conditional Distribution Parameters
#'
#' @param nm_result normal mixture model parameters
#' @param mean_residuals mean residuals by area
#' @param sigmae2est error variance
#' @param framework framework list with data and configuration
#'
#' @return conditional parameters by area
compute_conditional_params <- function(nm_result, mean_residuals, sigmae2est, framework) {
  if (!all(c("lambda", "mu", "sigma") %in% names(nm_result))) {
    stop("Invalid nm_result structure: Missing required components")
  }
  
  if (!nm_result$convergence) {
    warning("Using non-converged mixture model. Results may be unreliable.")
  }
  
  areas <- unique(as.character(framework$smp_data[, framework$smp_domains]))
  n_by_area <- table(as.character(framework$smp_data[, framework$smp_domains]))
  conditional_params <- vector("list", length(areas))
  names(conditional_params) <- areas
  
  for (area in areas) {
    n_i <- ifelse(area %in% names(n_by_area), n_by_area[area], 0)
    
    if (n_i == 0) {
      warning(sprintf("No observations for area: %s - Using defaults", area))
      conditional_params[[area]] <- list(
        weights = nm_result$lambda,
        means = nm_result$mu,
        variances = (nm_result$sigma)^2,
        n = 0,
        ebar = NA
      )
      next
    }
    
    ebar_i <- if (!is.na(mean_residuals[area])) {
      mean_residuals[area]
    } else {
      warning(sprintf("Missing residual for area: %s - Using mean", area))
      mean(mean_residuals, na.rm = TRUE)
    }
    
    if (any(is.na(nm_result$mu)) || any(is.na(nm_result$sigma))) {
      warning(sprintf("Invalid mixture params for area: %s - Using fallback", area))
      nm_result$mu <- rep(mean(mean_residuals, na.rm = TRUE), length(nm_result$lambda))
      nm_result$sigma <- rep(sd(mean_residuals, na.rm = TRUE), length(nm_result$lambda))
    }
    
    weights <- nm_result$lambda
    means <- nm_result$mu
    vars <- (nm_result$sigma)^2
    
    gamma <- vars / (vars + sigmae2est / n_i)
    cond_means <- gamma * ebar_i + (1 - gamma) * means
    cond_vars <- vars * sigmae2est / (vars + sigmae2est / n_i)
    
    cond_weights <- weights * dnorm(ebar_i, means, sqrt(vars + sigmae2est / n_i))
    cond_weights <- cond_weights / sum(cond_weights, na.rm = TRUE)
    
    conditional_params[[area]] <- list(
      weights = ifelse(is.na(cond_weights), weights, cond_weights),
      means = ifelse(is.na(cond_means), means, cond_means), 
      variances = ifelse(is.na(cond_vars), vars, cond_vars),
      n = n_i,
      ebar = ebar_i
    )
  }
  
  conditional_params
}

#' Calculate Indicators
#'
#' @param y response values
#' @param framework framework list with data and configuration
#' @param threshold poverty threshold
#'
#' @return list of indicators by domain
calculate_indicators <- function(y, framework, threshold) {
  domains <- if (!is.null(attributes(y)$domains)) {
    attributes(y)$domains
  } else {
    framework$pop_domains_vec[1:length(y)]
  }
  
  if (length(y) != length(domains)) {
    stop(sprintf("Length mismatch: values (%d) != domains (%d)", length(y), length(domains)))
  }
  
  y <- as.numeric(y)
  domain_values <- split(y, domains)
  domain_values <- lapply(domain_values, function(x) x[!is.na(x)])
  
  indicators <- list()
  
  indicators$Mean <- sapply(domain_values, function(x) {
    if (length(x) == 0) NA else mean(x, na.rm = TRUE)
  })
  
  # For threshold, apply function if it's a function, otherwise use as is
  indicators$Head_Count <- sapply(domain_values, function(x) {
    if (length(x) == 0) return(NA)
    # Calculate threshold - either from function or use provided value
    t <- if (is.function(threshold)) threshold(x) else threshold
    # Count values below threshold and divide by total count
    mean(x < t, na.rm = TRUE)
  })
  
  indicators$Poverty_Gap <- sapply(domain_values, function(x) {
    if (length(x) == 0) return(NA)
    t <- if (is.function(threshold)) threshold(x) else threshold
    gaps <- (t - x) / t
    gaps[gaps < 0] <- 0
    mean(gaps, na.rm = TRUE)
  })
  
  for (q in c(0.1, 0.25, 0.5, 0.75, 0.9)) {
    name <- if (q == 0.5) "Median" else paste0("Quantile_", q * 100)
    indicators[[name]] <- sapply(domain_values, function(x) {
      if (length(x) == 0) return(NA)
      quantile(x, probs = q, na.rm = TRUE)
    })
  }
  
  indicators$Gini <- sapply(domain_values, function(x) {
    if (length(x) < 2) return(NA)
    
    # Sort values in ascending order
    x <- sort(x, na.last = NA)
    n <- length(x)
    
    # Calculate Gini coefficient using the formula based on ordered values
    # G = (2 * sum(i*x_i) / (n * sum(x_i))) - (n+1)/n
    numerator <- sum(seq_along(x) * x)
    denominator <- sum(x)
    
    if (denominator <= 0) return(NA)  # Avoid division by zero
    
    g <- (2 * numerator) / (n * denominator) - (n + 1) / n
    
    # Ensure value is within [0,1] bounds
    return(max(0, min(1, g)))
  })
  
  indicators$Quintile_Share <- sapply(domain_values, function(x) {
    # Check if we have enough values for quintiles
    if (length(x) < 5) return(NA)
    
    # Remove negative and zero values to avoid issues
    x_pos <- x[x > 0]
    if (length(x_pos) < 5) return(NA)
    
    # Calculate quintiles
    quintiles <- quantile(x_pos, probs = c(0.2, 0.8), na.rm = TRUE)
    
    # Calculate mean income for top and bottom 20%
    top20 <- mean(x_pos[x_pos >= quintiles[2]], na.rm = TRUE)
    bottom20 <- mean(x_pos[x_pos <= quintiles[1]], na.rm = TRUE)
    
    # Avoid division by zero or negative ratio
    if (bottom20 <= 0) return(NA)
    
    ratio <- top20 / bottom20
    
    # Ensure ratio is positive and reasonable
    return(max(0, min(ratio, 100)))
  })
  
  # Handle custom indicators if provided
  if (!is.null(framework$custom_indicator)) {
    for (ind_name in names(framework$custom_indicator)) {
      indicators[[ind_name]] <- sapply(domain_values, function(x) {
        if (length(x) == 0) return(NA)
        t <- if (is.function(threshold)) threshold(x) else threshold
        # Calculate custom indicator with appropriate parameters
        tryCatch({
          framework$custom_indicator[[ind_name]](x, threshold = t)
        }, error = function(e) {
          # Try without threshold parameter if that fails
          tryCatch({
            framework$custom_indicator[[ind_name]](x)
          }, error = function(e2) {
            warning("Failed to calculate custom indicator: ", ind_name)
            NA
          })
        })
      })
    }
  }
  
  indicators
}

#' Back Transformation
#'
#' @param y transformed values
#' @param transformation transformation type
#' @param lambda transformation parameter
#' @param shift shift parameter
#' @param framework framework list with data and configuration
#' @param fixed formula for fixed effects
#'
#' @return back-transformed values
back_transformation <- function(y, transformation, lambda, shift, framework, fixed) {
  if (transformation == "box.cox") {
    if (abs(lambda) < 1e-6) {
      return(exp(y))
    } else {
      return((lambda * y + 1)^(1 / lambda))
    }
  } else if (transformation == "log") {
    return(exp(y))
  } else if (transformation == "no") {
    return(y)
  } else {
    warning("Unsupported transformation type: ", transformation, ", returning input values")
    return(y)
  }
}

#' Aggregate Weighted Mean
#'
#' @param df values
#' @param by grouping variable
#' @param w weights
#'
#' @return weighted means by group
aggregate_weighted_mean <- function(df, by, w) {
  agg <- stats::aggregate(
    x = df * w,
    by = by,
    FUN = sum
  )
  
  w_sums <- stats::aggregate(
    x = w,
    by = by,
    FUN = sum
  )
  
  agg[,2] <- agg[,2] / w_sums[,2]
  agg
}

#' Weighted Statistics
#'
#' @param x values
#' @param w weights
#'
#' @return list with sorted values, weights, and cumulative weights
weighted_stats <- function(x, w = NULL) {
  if(is.null(w)) w <- rep(1, length(x))
  complete <- complete.cases(x, w)
  x <- x[complete]
  w <- w[complete]
  w <- w / sum(w)
  ord <- order(x)
  x <- x[ord]
  w <- w[ord]
  cw <- cumsum(w)
  
  list(
    x = x,
    w = w,
    cw = cw
  )
}

#' Weighted Quantile
#'
#' @param x values
#' @param w weights
#' @param probs probabilities
#'
#' @return weighted quantiles
weighted_quantile <- function(x, w = NULL, probs = seq(0, 1, 0.25)) {
  stats <- weighted_stats(x, w)
  
  sapply(probs, function(p) {
    if(p <= 0) return(min(stats$x))
    if(p >= 1) return(max(stats$x))
    i <- max(which(stats$cw <= p))
    if(i < length(stats$x)) {
      x1 <- stats$x[i]
      x2 <- stats$x[i + 1]
      w1 <- stats$cw[i]
      w2 <- stats$cw[i + 1]
      x1 + (x2 - x1) * (p - w1) / (w2 - w1)
    } else {
      stats$x[i]
    }
  })
}

#' Compare NNEBP with standard EBP
#'
#' @param nnebp_obj NNEBP object
#' @param ebp_obj EBP object
#' @param indicators indicators to compare
#'
#' @return comparison data frame
compare_nnebp_ebp <- function(nnebp_obj, ebp_obj, indicators = "all") {
  # Validate inputs
  if (!inherits(nnebp_obj, "nnebp")) {
    stop("nnebp_obj must be an object of class 'nnebp'")
  }
  
  if (!inherits(ebp_obj, "ebp")) {
    stop("ebp_obj must be an object of class 'ebp'")
  }
  
  # Extract point estimates
  nnebp_ind <- nnebp_obj$ind
  ebp_ind <- ebp_obj$ind
  
  # Select indicators
  if (identical(indicators, "all")) {
    indicators <- intersect(names(nnebp_ind), names(ebp_ind))
  } else {
    indicators <- intersect(indicators, intersect(names(nnebp_ind), names(ebp_ind)))
  }
  
  if (length(indicators) == 0) {
    stop("No common indicators found between NNEBP and EBP results")
  }
  
  # Create comparison data frame
  result <- data.frame(Domain = names(nnebp_ind[[1]]))
  
  for (ind in indicators) {
    result[[paste0(ind, "_NNEBP")]] <- nnebp_ind[[ind]]
    result[[paste0(ind, "_EBP")]] <- ebp_ind[[ind]]
    
    # Add difference and relative difference
    result[[paste0(ind, "_Diff")]] <- nnebp_ind[[ind]] - ebp_ind[[ind]]
    result[[paste0(ind, "_RelDiff")]] <- (nnebp_ind[[ind]] - ebp_ind[[ind]]) / ebp_ind[[ind]] * 100
    
    # Add MSE if available
    if (!is.null(nnebp_obj$MSE) && !is.null(ebp_obj$MSE)) {
      result[[paste0(ind, "_MSE_NNEBP")]] <- nnebp_obj$MSE[[ind]]
      result[[paste0(ind, "_MSE_EBP")]] <- ebp_obj$MSE[[ind]]
      result[[paste0(ind, "_MSE_Diff")]] <- nnebp_obj$MSE[[ind]] - ebp_obj$MSE[[ind]]
    }
  }
  
  return(result)
}
