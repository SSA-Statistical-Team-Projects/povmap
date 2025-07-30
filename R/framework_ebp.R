# Internal documentation -------------------------------------------------------

# The function notation defines the notational framework for the EBP approach
# e.g. number of households in population or sample (per domain), distinction
# between in-sample and out-of-sample
# see Molina and Rao (2003) p.370-371


framework_ebp <- function(fixed, pop_data, pop_domains, pop_subdomains, smp_data, smp_domains,
                          smp_subdomains, threshold, custom_indicator = NULL, na.rm,
                          aggregate_to = NULL, weights, pop_weights, 
                          MSE_pop_weights, weights_type, benchmark_level, 
                          benchmark_weights,nlme_maxiter, nlme_tolerance, 
                          nlme_opt, nlme_optimmethod, nlme_method, nlme_mstol, 
                          nlme_returnobject, nlme_msmaxiter, rescale_weights,model_parameters,data.table,indicators) {

  # Reduction of number of variables
  mod_vars <- all.vars(fixed)
  mod_vars <- mod_vars[mod_vars != as.character(fixed[2])]

  if (!is.null(weights) && benchmark_weights==weights) {
    smp_vars <- c(as.character(fixed[2]), mod_vars, smp_domains, weights,
                  benchmark_level)
    if (!is.null(smp_subdomains)) {
      smp_vars <- c(smp_vars,smp_subdomains)
    }
  } else {
    smp_vars <- c(as.character(fixed[2]), mod_vars, smp_domains, weights,
                  benchmark_level, benchmark_weights)
    if (!is.null(smp_subdomains)) {
      smp_vars <- c(smp_vars,smp_subdomains)
    }
  }

  pop_vars <- c(mod_vars, pop_domains, pop_subdomains, aggregate_to, pop_weights,
                benchmark_level)
  if (!is.null(MSE_pop_weights)) {
  pop_vars <- c(pop_vars,MSE_pop_weights)
  }
  if (min(smp_vars %in% colnames(smp_data))==0) {
    stop(strwrap(prefix = " ", initial = "",
                paste0(smp_vars[smp_vars %in% colnames(smp_data)==F]))," are not in the sample dataset")
                  }
  smp_data <- smp_data[, smp_vars]
  weights <- weights
  pop_weights <- pop_weights
  fw_check1(
    pop_data = pop_data, mod_vars = mod_vars, pop_domains = pop_domains,
    pop_subdomains=pop_subdomains, smp_data = smp_data, aggregate_to = aggregate_to, fixed = fixed,
    smp_domains = smp_domains, smp_subdomains=smp_subdomains, threshold = threshold, weights = weights,
    pop_weights = pop_weights, benchmark_level = benchmark_level,
    benchmark_weights = benchmark_weights, weights_type = weights_type,
    rescale_weights = rescale_weights
  )

  if (min(pop_vars %in% colnames(pop_data))==0) {
    stop(strwrap(prefix = " ", initial = "",
                 paste0(pop_vars[pop_vars %in% colnames(pop_data)==F]))," are not in the population dataset")
  }
  pop_data <- pop_data[, pop_vars]
  # convert to dataframe if necessary 
  if ("tbl_df" %in% class(pop_data)) {
    pop_data <- as.data.frame(pop_data)
  }
  if ("tbl_df" %in% class(smp_data)) {
    smp_data <- as.data.frame(smp_data)
  }
  
  
  # Deletion of NA
  if (na.rm == TRUE) {
    pop_data <- na.omit(pop_data)
    smp_data <- na.omit(smp_data)
  } else if (any(is.na(pop_data)) || any(is.na(smp_data))) {
    stop(strwrap(prefix = " ", initial = "",
                 "EBP does not work with missing values. Set na.rm = TRUE in
                 function ebp."))
  }
  
  
  
  
  # rescale weights such that mean is equal to one within each domain 
  if (isTRUE(rescale_weights) && !is.null(weights)) {
    smp_data[,weights] <- smp_data[,weights] / ave(smp_data[,weights], smp_data[,smp_domains])
  }
  else if (!is.null(weights)) {
    smp_data[,weights] <- smp_data[,weights] / mean(smp_data[,weights])
  }
  
  # Order of domains
  pop_data <- pop_data[order(pop_data[[pop_domains]]), ]

  levels_tmp <- unique(pop_data[[pop_domains]])
  pop_data[[pop_domains]] <- factor(pop_data[[pop_domains]],
                                    levels = levels_tmp)
  pop_domains_vec <- pop_data[[pop_domains]]
  

  smp_data[[smp_domains]] <- factor(smp_data[[smp_domains]],
                                    levels = levels_tmp)

 
  #This keeps common subdomains in sample, which is causing several NAs.  
 if (!is.null(pop_subdomains)) {
pop_data <- pop_data[order(pop_data[[pop_subdomains]]), ]
levels_subdom_tmp <- unique(pop_data[[pop_subdomains]])
pop_data[[pop_subdomains]] <- factor(pop_data[[pop_subdomains]],
                                    levels = levels_subdom_tmp)
  pop_subdomains_vec <- pop_data[[pop_subdomains]]
  # levels option in smp_data is intentionally omitted to not restrict subdomains to set of population values 
  # If we only take observations with matching population agebs it will create  
  # missing values in the subdomain id variable which will mess up model estimation 
  smp_data[[smp_subdomains]] <- factor(smp_data[[smp_subdomains]])
  }
  
  
  
  
  
  if (is.null(aggregate_to)) {
    aggregate_to_vec <- NULL
  } else {
    levels_tmp <- unique(pop_data[[aggregate_to]])
    pop_data[[aggregate_to]] <- factor(pop_data[[aggregate_to]],
                                       levels = levels_tmp)
    aggregate_to_vec <- pop_data[[aggregate_to]]
  }

  rm(levels_tmp)
  smp_data <- smp_data[order(smp_data[[smp_domains]]), ]


  smp_domains_vec <- smp_data[[smp_domains]]
  smp_domains_vec <- droplevels(smp_domains_vec)
  
  smp_subdomains_vec <- NULL 
  pop_subdomains_vec <- NULL 
  both_subdomains_vec <- NULL 
  if (!is.null(smp_subdomains) && !is.null(pop_subdomains)) {
    pop_subdomains_vec <- pop_data[[pop_subdomains]]
    smp_subdomains_vec <- smp_data[[smp_subdomains]]
    both_subdomains_vec <- smp_subdomains_vec[smp_subdomains_vec %in% pop_subdomains_vec]
  }
  
  
  

  fw_check2(
    pop_domains = pop_domains, pop_domains_vec = pop_domains_vec,
    smp_domains = smp_domains, smp_domains_vec = smp_domains_vec,
    aggregate_to = aggregate_to, aggregate_to_vec = aggregate_to_vec
  )


  # Number of households in population
  N_pop <- length(pop_domains_vec)
  # Number of households in sample
  N_smp <- length(smp_domains_vec)
  # Number of out-of-sample households
  N_unobs <- N_pop - N_smp
  # Number of domains in the population
  N_dom_pop <- length(unique(pop_domains_vec))
  # Number of subdomains in the population 
  N_subdom_pop <- length(unique(pop_subdomains_vec))
  # Number of domains in the population on aggregated level
  N_dom_pop_agg <- length(unique(aggregate_to_vec))
  # Number of domains in the sample
  N_dom_smp <- length(unique(smp_domains_vec))
  # Number of out-of-sample domains
  N_dom_unobs <- N_dom_pop - N_dom_smp
  # Number of subdomains in sample 
  N_subdom_smp <- length(unique(both_subdomains_vec))
  # Number of out-of-sample subdomains
  N_subdom_unobs <- N_subdom_pop - N_subdom_smp
  # Number of households in population per domain
  n_pop <- as.vector(table(pop_domains_vec))
  # NUmber of households in population per subdomain 
  n_pop_subdom <- as.vector(table(pop_subdomains_vec))
  # Number of households in sample per domain
  smp_domains_vec_tmp <- as.numeric(smp_domains_vec)
  n_smp <- as.vector(table(smp_domains_vec_tmp))
  smp_subdomains_vec_tmp <- as.numeric(smp_subdomains_vec)
  n_smp_subdom <- as.vector(table(smp_subdomains_vec_tmp))
  
  
  
  # Indicator variables that indicate if domain is in- or out-of-sample
  obs_dom <- pop_domains_vec %in% unique(smp_domains_vec)
  dist_obs_dom <- unique(pop_domains_vec) %in% unique(smp_domains_vec)
  obs_subdom <- pop_subdomains_vec %in% unique(smp_subdomains_vec)
  dist_obs_subdom <- unique(pop_subdomains_vec) %in% unique(smp_subdomains_vec)
  
  obs_smp_dom <- smp_domains_vec %in% unique(pop_domains_vec)
  dist_obs_smp_dom <- unique(smp_domains_vec) %in% unique(pop_domains_vec)
  
  obs_smp_subdom <- smp_subdomains_vec %in% unique(pop_subdomains_vec)
  dist_obs_smp_subdom <- unique(smp_subdomains_vec) %in% unique(pop_subdomains_vec)
  
  fw_check3(
    obs_dom = obs_dom, dist_obs_dom = dist_obs_dom, pop_domains = pop_domains,
    smp_domains = smp_domains
  )

  indicator_list <- list(
    Mean = function(y, pop_weights, threshold) {
      t(weighted.mean(y, pop_weights))
    },
    Head_Count = function(y, pop_weights, threshold) {
       t(weighted.mean(y < threshold, pop_weights))
    },
    Poverty_Gap = function(y, pop_weights, threshold) {
      sum((1 - (y[y < threshold]/threshold)) * pop_weights[y < threshold])/
        sum(pop_weights)
    },
    Poverty_Severity = function(y, pop_weights, threshold) {
      sum((1 - y[y < threshold]/threshold)^2 * pop_weights[y < threshold])/
        sum(pop_weights)
    },
    Gini = function(y, pop_weights, threshold) {
        n <- length(y)
        pop_weights <- pop_weights[order(y)]
        y <- sort(y)
        auc <- sum((cumsum(c(0, (y * pop_weights)[1:(n-1)])) +
                      ((y * pop_weights) / 2)) * pop_weights)
        auc <- (auc / sum(pop_weights)) / sum((y * pop_weights))
        G <- 1 - (2* auc)
        return(G)
    },
    Quintile_Share = function(y, pop_weights, threshold) {
      quant14 <- wtd.quantile(x = y, weights = pop_weights, probs = c(0.2, 0.8))

      iq1 <- y <= quant14[1]
      iq4 <- y > quant14[2]
      t((sum(pop_weights[iq4] * y[iq4]) / sum(pop_weights[iq4])) /
           (sum(pop_weights[iq1] * y[iq1]) / sum(pop_weights[iq1])))
    },
    Quantiles = function(y, pop_weights, threshold) {
      if(length(unique(pop_weights)) == 1 & 1 %in% unique(pop_weights)){
        t(quantile(x = y, probs = c(0, .10, .25, .5, .75, .9, 1)))
      }else{
        t(wtd.quantile(x = y, weights = pop_weights,
                       probs = c(0,.10, .25, .5, .75, .9, 1)))
      }
    }
  )

  if (data.table==TRUE) {
    indicator_list[["Mean"]] <-   
          function(y,pop_weights,threshold=NULL) {
      weighted.mean(x=y,w=pop_weights)
    }
    
    indicator_list[["Head_Count"]] <- 
     function (y,threshold,pop_weights) {
      # write functions to operate on vectors and then use data.table subsetting 
      weighted.mean(y<threshold,w=pop_weights)
     }
    
    indicator_list[["Quintile_Share"]] <- 
     function(y, pop_weights, threshold) {
        quant14 <- wtd.quantile(x = y, weights = pop_weights, probs = c(0.2, 0.8))
        
        iq1 <- y <= quant14[1]
        iq4 <- y > quant14[2]
        (sum(pop_weights[iq4] * y[iq4]) / sum(pop_weights[iq4])) /
            (sum(pop_weights[iq1] * y[iq1]) / sum(pop_weights[iq1]))
      }
    
    indicator_list[["Quantiles"]] <-
      function(y, pop_weights, threshold) {
        if(length(unique(pop_weights)) == 1 & 1 %in% unique(pop_weights)){
          quantile(x = y, probs = c(.10, .25, .5, .75, .9))
        }else{
          wtd.quantile(x = y, weights = pop_weights,
                         probs = c(.10, .25, .5, .75, .9))
        }
      }
      
    
  } # close data.table functions 
  
  
  indicator_names <- c(
    "Mean",
    "Head_Count",
    "Poverty_Gap",
    "Poverty_Severity",
    "Gini",
    "Quintile_Share",
    "Quantile_10",
    "Quantile_25",
    "Median",
    "Quantile_75",
    "Quantile_90"
  )

  function_names <- c(
    "Mean",
    "Head_Count",
    "Poverty_Gap",
    "Poverty_Severity",
    "Gini",
    "Quintile_Share",
    "Quantiles"
  )
  

  
  if (!is.null(indicators)) {
    keepthese <- which(function_names %in% indicators)
    indicator_list <- indicator_list[keepthese]
    if (7 %in% keepthese) {
      keepthese <- c(keepthese,8,9,10,11)
    }
    indicator_names <- indicator_names[keepthese]
  }

  if (!is.null(custom_indicator) && length(custom_indicator) > 0) {
    for(i in 1:length(custom_indicator)) {
      formals(custom_indicator[[i]]) <- alist(y=, pop_weights=, threshold=)
    }

    indicator_list <- c(indicator_list, custom_indicator)
    indicator_names <- c(indicator_names, names(custom_indicator))
  }

  if (is.null(threshold)) {
    threshold <- 0.6 * median(smp_data[[paste(fixed[2])]])
    message(strwrap(prefix = " ", initial = "",
                    paste0("The threshold for the HCR and the PG is
                          automatically set to 60% of the median of the
                          dependent variable and equals ", threshold)))
  }

  return(list(
    pop_data = pop_data,
    pop_domains_vec = pop_domains_vec,
    pop_domains = pop_domains, 
    pop_subdomains = pop_subdomains, 
    pop_subdomains_vec = pop_subdomains_vec, 
    smp_data = smp_data,
    smp_domains_vec = smp_domains_vec,
    smp_domains = smp_domains,
    smp_subdomains = smp_subdomains,
    smp_subdomains_vec = smp_subdomains_vec, 
    both_subdomains_vec = both_subdomains_vec, 
    #smp_subdomains_vec_both = smp_subdomains_vec_both, 
    aggregate_to = aggregate_to,
    aggregate_to_vec = aggregate_to_vec,
    N_pop = N_pop,
    N_smp = N_smp,
    N_unobs = N_unobs,
    N_dom_pop = N_dom_pop,
    N_dom_pop_agg = N_dom_pop_agg,
    N_subdom_pop = N_subdom_pop, 
    N_dom_smp = N_dom_smp,
    N_dom_unobs = N_dom_unobs,
    N_subdom_smp = N_subdom_smp,
    N_subdom_unobs = N_subdom_unobs,
    n_pop = n_pop,
    n_smp = n_smp,
    n_pop_subdom = n_pop_subdom, 
    n_smp_subdom = n_smp_subdom,
    obs_dom = obs_dom,
    obs_subdom = obs_subdom, 
    dist_obs_dom = dist_obs_dom,
    dist_obs_subdom = dist_obs_subdom,
    obs_smp_dom = obs_smp_dom,
    dist_obs_smp_dom = dist_obs_smp_dom,
    obs_smp_subdom = obs_smp_subdom,
    dist_obs_smp_subdom = dist_obs_smp_subdom,
    indicator_list = indicator_list,
    indicator_names = indicator_names,
    threshold = threshold,
    weights = weights,
    benchmark_weights = benchmark_weights,
    pop_weights = pop_weights,
    MSE_pop_weights = MSE_pop_weights, 
    weights_type = weights_type,
    nlme_maxiter = nlme_maxiter,
    nlme_tolerance = nlme_tolerance,
    nlme_opt = nlme_opt, 
    nlme_optimmethod = nlme_optimmethod, 
    nlme_msmaxiter = nlme_msmaxiter, 
    nlme_mstol = nlme_mstol, 
    nlme_returnobject = nlme_returnobject, 
    nlme_method = nlme_method,
    model_parameters = model_parameters,
    data.table=data.table,
    indicators=indicators
  ))
}

# A simple scaling function for weights
scaler <- function(x, region){
  average_x <- tapply(x, INDEX = region, FUN = mean, na.rm = TRUE)
  df <- merge(data.frame(x, region), data.frame(names(average_x), average_x),
              by.x = "region", by.y = "names.average_x.")
  return(df$x / df$average_x)
}
