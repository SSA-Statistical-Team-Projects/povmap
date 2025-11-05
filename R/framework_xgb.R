framework_xgb <-function(fixed,
                         smp_data,
                         pop_data,
                         smp_weights = NULL,
                         pop_weights = NULL,
                         domains,
                         transformation,
                         conf_level,
                         sub_domains,
                         benchmark_level, 
                         benchmark_weights, 
                         na.rm) {

  # Data preparation
  # Splitting the fixed string to extract outcome and covariates
  split <- strsplit(as.character(fixed[-1]), "~", fixed = TRUE)
  outcome <- trimws(split[[1]][1])
  covariates <- trimws(strsplit(trimws(split[[2]]), "\\+")[[1]])

  #checks 
  if (!all(covariates %in% colnames(pop_data))) {
    missing_vars <- covariates[(!covariates %in% colnames(pop_data))]
    stop(paste("Variables",missing_vars,"not present in population dataframe"))
  }
  if (!domains %in% colnames(pop_data)) {
    stop(paste("Domain identifier",domains,"not present in population dataframe"))
  }
  if (!sub_domains %in% colnames(pop_data)) {
    stop(paste("Subdomain identifier",sub_domains,"not present in population dataframe"))
  }
  if (!outcome %in% colnames(smp_data)) {
    stop(paste("Outcome",outcome,"not present in sample dataframe"))
  }
  

  # Deletion of NA
  if (na.rm == TRUE) {
    pop_data <- na.omit(pop_data)
    smp_data <- na.omit(smp_data)
  } else if (any(is.na(pop_data)) || any(is.na(smp_data))) {
    stop(strwrap(prefix = " ", initial = "",
                 "XGB does not work with missing values. Set na.rm = TRUE in
                 function xgb."))
  }

  # Extracting relevant subsets of data
  X_smp <- smp_data[, c(covariates,domains,sub_domains)]
  
  if (!is.null(smp_weights) && benchmark_weights==smp_weights) {
     X_smp <- cbind(X_smp,smp_data[,benchmark_level])
  }
  else {
     X_smp <- cbind(X_smp,smp_data[,c(benchmark_level,benchmark_weights)])
  }
    
  
  
  
  
  Y_smp <- smp_data[, outcome]
  #X_pop <- pop_data[, c(covariates,domains,sub_domains,pop_weights)]
  X_pop <- pop_data[, c(covariates,domains,sub_domains)]
  # add pop weights if they are not already in thecovariates 
  if (!pop_weights %in% covariates) {
    X_pop <- data.frame(X_pop,pop_weights)
  }
  
  # Handling sample and population weights
    smp_weights_name <- smp_weights 
    if (!is.null(smp_weights)) {
    smp_weights <- smp_data[, smp_weights]
  } else {
    smp_weights <- rep(1, length = length(Y_smp))
  }

    pop_weights_name <- pop_weights 
  if (!is.null(pop_weights)) {
    pop_weights <- pop_data[, pop_weights]
  } else {
    pop_weights <- rep(1, length = nrow(pop_data))
  }


  

  # Determining domains in sample and population
  in_smp <- unique(smp_data[[domains]])
  total_dom <- unique(pop_data[[domains]])
  out_smp <- !total_dom %in% in_smp

  # Storing summary information
    N_smp = length(smp_data[[domains]])
    N_pop = length(pop_data[[domains]])
    domains_out = sum(out_smp)
    domains_in = length(in_smp)
    domains_total = length(total_dom)
    ni_smp = table(smp_data[[domains]])
    ni_pop = table(pop_data[[domains]])
    domains = domains
    sub_domains = sub_domains
    outcome = outcome
    smp_weights_name = smp_weights_name 
    pop_domains_vec <- pop_data[[domains]]

  # Check
  xgb_check1(
    transformation = transformation,
    Y_smp = Y_smp,
    X_smp = X_smp,
    X_pop = X_pop,
    smp_weights = smp_weights,
    pop_weights = pop_weights,
    conf_level = conf_level,
    domains = domains,
    sub_domains = sub_domains
  )

  # Transformation
  # Applying specified transformations to Y_smp
  if (!is.null(transformation)) {
    if (transformation == "arcsin") {
      Y_smp <- asin(sqrt(Y_smp))
    } else if (transformation == "log") {
      Y_smp <- log(Y_smp)
    }
  }

  return(list(Y_smp = Y_smp,
              X_smp = X_smp,
              X_pop = X_pop,
              smp_data = smp_data, 
              pop_data = pop_data, 
              smp_weights_vec = smp_weights,
              pop_weights_vec = pop_weights,
              N_smp = N_smp, 
              N_pop = N_pop, 
              domains_out = domains_out, 
              domains_in = domains_in, 
              domains_total = domains_total, 
              ni_smp = ni_smp, 
              ni_pop = ni_pop, 
              domains = domains, 
              sub_domains = sub_domains,
              outcome=outcome, 
              smp_weights = smp_weights_name,
              pop_weights = pop_weights_name, 
              covariates = covariates,
              benchmark_weights = benchmark_weights, 
              smp_domains = domains, 
              pop_domains_vec = pop_domains_vec 
              ))
}