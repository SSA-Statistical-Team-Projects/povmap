framework_xgb <-function(fixed,
                         smp_data,
                         pop_data,
                         smp_weights = NULL,
                         pop_weights = NULL,
                         domains,
                         transformation,
                         conf_level,
                         sub_domains = NULL,
                         na.rm) {

  # Data preparation
  # Splitting the fixed string to extract outcome and covariates
  split <- strsplit(as.character(fixed[-1]), "~", fixed = TRUE)
  outcome <- trimws(split[[1]][1])
  covariates <- trimws(strsplit(trimws(split[[2]]), "\\+")[[1]])



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
  X_smp <- smp_data[, c(domains,subdomains,covariates)]
  Y_smp <- smp_data[, c(domains,outcome)]
  
  # Handling sample and population weights
  if (!is.null(smp_weights)) {
    smp_weights <- smp_data[, smp_weights]
  } else {
    smp_weights <- rep(1, length = nrow(Y_smp))
  }

  if (!is.null(pop_weights)) {
    pop_weights <- pop_data[, pop_weights]
  } else {
    pop_weights <- rep(1, length = nrow(pop_data))
  }

  X_pop <- pop_data[, c(domains,covariates)]

  # Determining domains in sample and population
  in_smp <- unique(smp_data[[domains]])
  total_dom <- unique(pop_data[[domains]])
  out_smp <- !total_dom %in% in_smp

  # Storing summary information
  saeinfo <- list(
    N_smp = length(smp_data[[domains]]),
    N_pop = length(pop_data[[domains]]),
    domains_out = sum(out_smp),
    domains_in = length(in_smp),
    domains_total = length(total_dom),
    ni_smp = table(smp_data[[domains]]),
    ni_pop = table(pop_data[[domains]]),
    domains = domains
  )

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
      Y_smp[,2] <- asin(sqrt(Y_smp[,2]))
    } else if (transformation == "log") {
      Y_smp[,2] <- log(Y_smp[,2])
    }
  }

  return(list(saeinfo = saeinfo,
              Y_smp = Y_smp,
              X_smp = X_smp,
              X_pop = X_pop,
              smp_weights = smp_weights,
              pop_weights = pop_weights))
}
