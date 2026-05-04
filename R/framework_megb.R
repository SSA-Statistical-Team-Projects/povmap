framework_megb <- function(fixed,
                           smp_data,
                           pop_data,
                           smp_weights = NULL,
                           pop_weights = NULL,
                           domains,
                           transformation,
                           conf_level,
                           na.rm,
                           benchmark_weights = NULL) {

  # Parse formula
  split      <- strsplit(as.character(fixed[-1]), "~", fixed = TRUE)
  outcome    <- trimws(split[[1]][1])
  cov_all    <- trimws(strsplit(trimws(split[[2]]), "\\+")[[1]])

  # Domain is the LMM grouping variable, not a GB predictor
  covariates <- setdiff(cov_all, domains)

  # Validation
  if (!outcome %in% colnames(smp_data))
    stop(paste("Outcome", outcome, "not present in sample dataframe"))
  if (!domains %in% colnames(smp_data))
    stop(paste("Domain identifier", domains, "not present in sample dataframe"))
  if (!domains %in% colnames(pop_data))
    stop(paste("Domain identifier", domains, "not present in population dataframe"))
  missing_cov <- setdiff(covariates, colnames(pop_data))
  if (length(missing_cov) > 0)
    stop(paste("Variables not present in population dataframe:",
               paste(missing_cov, collapse = ", ")))

  valid_transforms <- c("no", "log", "log.shift", "logistic", "arcsin", "poisson")
  if (!transformation %in% valid_transforms)
    stop(paste("transformation must be one of:",
               paste(valid_transforms, collapse = ", ")))
  if (!is.numeric(conf_level) || conf_level <= 0 || conf_level >= 1)
    stop("conf_level must be a number strictly between 0 and 1")

  # NA handling
  if (na.rm) {
    smp_data <- na.omit(smp_data)
    pop_data <- na.omit(pop_data)
  } else {
    if (any(is.na(smp_data[[domains]])) || any(is.na(pop_data[[domains]])))
      stop(paste("Domain variable", domains,
                 "contains missing values. Set na.rm = TRUE to remove them."))
  }

  # Filter zero-weight population observations
  pop_weights_name <- pop_weights
  if (!is.null(pop_weights)) {
    pop_data         <- pop_data[pop_data[[pop_weights]] > 0, ]
    pop_weights_vec  <- pop_data[[pop_weights]]
  } else {
    pop_weights_vec  <- rep(1, nrow(pop_data))
  }

  smp_weights_name <- smp_weights
  if (!is.null(smp_weights)) {
    smp_weights_vec <- smp_data[[smp_weights]]
  } else {
    smp_weights_vec <- rep(1, nrow(smp_data))
  }

  # Domain summary info
  in_smp    <- unique(smp_data[[domains]])
  total_dom <- unique(pop_data[[domains]])

  ni_pop_tbl <- table(pop_data[[domains]])

  list(
    Y_smp             = smp_data[[outcome]],
    smp_data          = smp_data,
    pop_data          = pop_data,
    smp_weights_vec   = smp_weights_vec,
    pop_weights_vec   = pop_weights_vec,
    N_smp             = nrow(smp_data),
    N_pop             = nrow(pop_data),
    domains_out       = sum(!total_dom %in% in_smp),
    domains_in        = length(in_smp),
    domains_total     = length(total_dom),
    ni_smp            = table(smp_data[[domains]]),
    ni_pop            = ni_pop_tbl,
    n_pop             = as.vector(ni_pop_tbl),
    domains           = domains,
    smp_domains       = domains,
    outcome           = outcome,
    covariates        = covariates,
    smp_weights       = smp_weights_name,
    pop_weights       = pop_weights_name,
    benchmark_weights = if (!is.null(benchmark_weights)) benchmark_weights
                        else smp_weights_name,
    aggregate_to      = NULL,
    pop_domains_vec   = pop_data[[domains]]
  )
}
