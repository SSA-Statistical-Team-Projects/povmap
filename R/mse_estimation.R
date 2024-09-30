# Internal documentation -------------------------------------------------------

# MSE estimation - parametric bootstrap procedure

# Function parametric_bootstrap conducts the MSE estimation defined in function
# mse_estim (see below)
# The parametric boostrap approach can be find in Molina and Rao (2010) p. 376

parametric_bootstrap <- function(framework,
                                 point_estim,
                                 fixed,
                                 transformation,
                                 interval = c(-1, 2),
                                 L,
                                 B,
                                 boot_type,
                                 parallel_mode,
                                 cpus,
                                 benchmark,
                                 benchmark_type,
                                 benchmark_level) {


  message("\r", "Bootstrap started                                            ")
  if (boot_type == "wild") {
    res_s <- residuals(point_estim$model)
    fitted_s <- fitted(point_estim$model, level = 1)
  } else {
    res_s <- NULL
    fitted_s <- NULL
  }
  start_time <- Sys.time()
  if (cpus > 1) {
    cpus <- min(cpus, parallel::detectCores())
    parallelMap::parallelStart(
      mode = parallel_mode,
      cpus = cpus, show.info = F, logging = F
    )

    if (parallel_mode == "socket") {
      parallel::clusterSetRNGStream()
    }
    
    # check to see if current library matches slave library locations 
    # if not there will be a bug unless we alter the slave library locations 
    slavelibs <- parallel::clusterEvalQ(, .libPaths())
    if (!max(.libPaths() %in% slavelibs[[1]])) {
      libname <- .libPaths()
      parallelMap:::parallelExport("libname")
      parallel:::clusterEvalQ(, .libPaths(libname))
    }
    
    parallelMap::parallelLibrary("nlme")
    
    mses <- simplify2array(parallelMap::parallelLapply(
      xs              = seq_len(B),
      fun             = mse_estim_wrapper,
      B               = B,
      framework       = framework,
      lambda          = point_estim$optimal_lambda,
      shift           = point_estim$shift_par,
      model_par       = point_estim$model_par,
      gen_model       = point_estim$gen_model,
      fixed           = fixed,
      transformation  = transformation,
      interval        = interval,
      L               = L,
      res_s           = res_s,
      fitted_s        = fitted_s,
      start_time      = start_time,
      boot_type       = boot_type,
      benchmark       = benchmark,
      benchmark_type  = benchmark_type,
      benchmark_level = benchmark_level
    ))
    parallelMap::parallelStop()
  } else {

    mses <- simplify2array(lapply(
      X = seq_len(B),
      FUN = mse_estim_wrapper,
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
      boot_type = boot_type,
      benchmark = benchmark,
      benchmark_type = benchmark_type,
      benchmark_level = benchmark_level
    ))
  }

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




# mse_estim (only internal) ----------------------------------------------------

# The mse_estim function defines all parameters and estimations which have to
# be replicated B times for the Parametric Bootstrap Approach.
# See Molina and Rao (2010) p. 376

mse_estim <- function(framework,
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
                      boot_type,
                      benchmark,
                      benchmark_type,
                      benchmark_level) {



  # The function superpopulation returns an income vector and a temporary
  # variable that passes the random effect to generating bootstrap populations
  # in bootstrap_par.

  if (is.null(framework$MSE_pop_weights)) { 
  if (boot_type == "wild") {
    superpop <- superpopulation_wild(
      framework = framework,
      model_par = model_par,
      gen_model = gen_model,
      lambda = lambda,
      shift = shift,
      transformation = transformation,
      res_s = res_s,
      fitted_s = fitted_s,
      fixed = fixed
    )
  } else if (is.null(framework$smp_subdomains) && is.null(framework$pop_subdomains)) {
    superpop <- superpopulation(
      framework = framework,
      model_par = model_par,
      gen_model = gen_model,
      lambda = lambda,
      shift = shift,
      transformation = transformation,
      fixed = fixed
    )
  } else { # two fold model 
    superpop <- superpopulation_2f(
      framework = framework,
      model_par = model_par,
      gen_model = gen_model,
      lambda = lambda,
      shift = shift,
      transformation = transformation,
      fixed = fixed
    )
  }
    
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
  vu_tmp <- superpop$vu_tmp 
  eta_tmp <- superpop$eta_tmp 
  
} # close if no MSE pop weighting 
  
  else {
    
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
    true_indicators_weighted <- true_indicators_weighted(
      framework = framework,
      model_par = model_par, 
      gen_model = gen_model,
      lambda = lambda, 
      shift = shift,
      transformation = transformation,
      fixed = fixed
    )
    true_indicators <- true_indicators_weighted$true_indicators
    vu_tmp<- true_indicators_weighted$vu_tmp
  }
  
  
  
  
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
        add_bench <- as.data.frame(true_indicators[, names(benchmark)[-1]])
        if (!is.null(dim(add_bench))) {
          colnames(add_bench) <- c(paste0(names(benchmark)[-1],"_bench"))
        } else {
          names(add_bench) <- c(paste0(names(benchmark)[-1],"_bench"))
        }
      }
    }
    true_indicators <- as.matrix(cbind(true_indicators, add_bench))
  }

  # The function bootstrap_par returns a sample that can be given into the
  # point estimation to get predictors of the indicators that can be compared
  # to the "truth".

  if (boot_type == "wild") {
    bootstrap_sample <- bootstrap_par_wild(
      fixed = fixed,
      transformation = transformation,
      framework = framework,
      model_par = model_par,
      lambda = lambda,
      shift = shift,
      vu_tmp = vu_tmp,
      res_s = res_s,
      fitted_s = fitted_s
    )
  } else {
    bootstrap_sample <- bootstrap_par(
      fixed = fixed,
      transformation = transformation,
      framework = framework,
      model_par = model_par,
      lambda = lambda,
      shift = shift,
      vu_tmp = vu_tmp
    )
  }

  framework$smp_data <- bootstrap_sample

  # Prediction of indicators with bootstap sample.
  bootstrap_point_estim <- as.matrix(point_estim(
    fixed = fixed,
    transformation =
      transformation,
    interval = interval,
    L = L,
    framework = framework
  )[[1]][, -1])

  colnames(bootstrap_point_estim) <- framework$indicator_names
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

# True_indicators_weighted function ---------------------------------------------
#This returns "true" superpopulation Mean and Headcount estimated from the population data when MSE_pop_Weight=TRUE 
true_indicators_weighted <- function(framework,model_par,gen_model,lambda,shift,transformation,fixed) {
  
if (is.null(framework$smp_subdomains) && is.null(framework$pop_subdomains)) {  # one fold model 
  eta_pop <- rep(0,framework$N_pop)
  framework$obs_subdom <- rep(0,framework$N_pop)
} 
  else { # 2 fold model 
    #sigmau2est <- model_par$sigmau2est 
    #sigmae2est <- model_par$sigmae2est
    #sigmah2est <- model_par$sigmah2est 
    eta_tmp <- rnorm(framework$N_subdom_pop,0,sqrt(model_par$sigmah2est)) # will be zero for one fold model 
    eta_pop <- rep(eta_tmp, framework$n_pop_subdom)
  }
  
  
  
  #if (framework$MSE_random_variance==TRUE) {
    # Even though the two error terms are nindependent, the estimaes of their variances are correlated 
    # We define X ~N(MuX,s2X) and Y=A(X-MuX)+B, and then Y ~ N(b,A^2*s2X+s2B) 
    # so s2b should be equal to s2Y-A^2*s2x 
    # we have an estimate of Cov(lnsigmau2est)=model_par$cov_sigma2est = A 
    # If we define lnsigmae2est as Y and lnsigmaeu2est as X, then setting  s2B = s2Y-s2x*A2 will generate s2Y=A^2S2X+S2Y-S2X*A2 
    #lnsigmau2est<-rnorm(n=1,mean=log(sigmau2est),sd=sqrt(model_par$var_lnsigmau2est))
    #lnsigmae2est <- (lnsigmau2est-log(sigmae2est)*model_par$cov_sigma2est+rnorm(n=1,mean=log(sigmae2est),sd=sqrt(model_par$var_lnsigmae2est-model_par$var_lnsigmau2est*model_par$cov_sigma2est^2)))
    #sigmau2est <- exp(lnsigmau2est)
    #sigmae2est <- exp(lnsigmae2est)
    # first figure out variance/covariance matrix using delta method 
    #var_sigmau2 <- sigmau2est^2*model_par$var_lnsigmau2est  # = exp(lnsigmau2est)^2*var(lnsigmau2est)
    #var_sigmae2 <- sigmae2est^2*model_par$var_lnsigmae2est  # = exp(lnsigmau2est)^2*var(lnsigmau2est)
    #cov_sigmau2_sigmae2 <- sigmau2est*sigmae2est*model_par$cov_sigma2est
    # now add on noise 
    #sigmau2est <- rnorm(n=1,mean=sigmau2est,sd=sqrt(var_sigmau2))
    #sigmae2est  <- (sigmau2est-model_par$sigmau2est)*cov_sigmau2_sigmae2+rnorm(n=1,mean=sigmae2est,sd=sqrt(var_sigmae2-var_sigmau2*cov_sigmau2_sigmae2^2))
    
    # method 3 - this doesn't change anything 
      #sigmae2est <- exp(log(sigmae2est)+0.5*model_par$var_lnsigmae2est)
      #sigmau2est <- exp(log(sigmau2est)+0.5*model_par$var_lnsigmau2est)
    
    # method 4 = simple version of method 1, with no covariance  
    #lnsigmau2est<-rnorm(n=1,mean=log(sigmau2est),sd=sqrt(model_par$var_lnsigmau2est))
    #lnsigmae2est<-rnorm(n=1,mean=log(sigmae2est),sd=sqrt(model_par$var_lnsigmae2est))
    #sigmau2est <- exp(lnsigmau2est)
    #sigmae2est <- exp(lnsigmae2est)
  #}
  
  
  # draw new superpopulation random effect
  

    vu_tmp <- rnorm(framework$N_dom_pop, 0, sqrt(model_par$sigmau2est))
    vu_pop <- rep(vu_tmp, framework$n_pop)
   
    

 
  
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
    pop_weights_vec <- rep(1, framework$N_pop)
  }
  
  
  #  Draw epsilon 
  var_eps <- vector(length=framework$N_pop)
  var_eps[(framework$obs_dom & framework$obs_subdom)] <- (model_par$sigmae2est) # observed subdomain within observed domain 
  var_eps[(framework$obs_dom & !framework$obs_subdom)] <- (model_par$sigmae2est+model_par$sigmah2est) # unobserved subdomain within observed subdomain 
  var_eps[!framework$obs_dom] <- (model_par$sigmae2est+model_par$sigmah2est+model_par$sigmau2est) # unobserved domain 
  
  eps <- rnorm(framework$N_pop, 0, sqrt(var_eps))

  
  true_indicators_weighted<-matrix(nrow = N_dom_pop_tmp, 
                                   ncol = length(framework$indicator_names))
  colnames(true_indicators_weighted) <- framework$indicator_names 
 #Do Mean calculation 
#1. back-transform draw   
  Y_pop_b <- gen_model$mu_fixed + vu_pop + eta_pop + eps
  Y_pop_b <- back_transformation(y=Y_pop_b,transformation=transformation,lambda=lambda,shift=shift,framework=framework)
#2. Calulate expected value of transformed XB+mu
  Y_pop_mu <- gen_model$mu_fixed + vu_pop + eta_pop 
  if ("Mean" %in% framework$indicator_names) {
    EY_pop_mu <- expected_untransformed_mean(Y_pop_mu,var=var_eps,transformation=transformation,lambda=lambda)
#3. Scale down implied residual to simulate taking draws over repeated observations   
  Y_pop_b <- EY_pop_mu+((Y_pop_b-Y_pop_mu)/sqrt(framework$pop_data[,framework$MSE_pop_weights]))
    true_indicators_weighted[,"Mean"] <- mapply(FUN=weighted.mean, x=split(Y_pop_b, pop_domains_vec_tmp),w=split(pop_weights_vec,pop_domains_vec_tmp))
    # Note that for no transformation case, true_indicators Y_pop_b ~ N(XB+mu, var(epsilon)/MSE_pop_weights) 
    # for log transformation case, draws log normal, calculates the implied residual in levels, and scales that down 
  }
  
  if ("Head_Count" %in% framework$indicator_names) {
# Do headcount calculation for population
# 1. Find poverty probability of each obserbation
    p_pov <- vector(length = framework$N_pop)
    p_pov <- expected_head_count(mu=Y_pop_mu,threshold=framework$threshold,var=var_eps,transformation=transformation,lambda=lambda,shift=shift)
# draw from binomial distribution, then divide by total number of trials      
    pov <- mapply(FUN=rbinom,n=1,size=framework$pop_data[,framework$MSE_pop_weights],prob=p_pov)
    pov <- as.vector(pov)/framework$pop_data[,framework$MSE_pop_weights] #divide by total number of trials 
# take weighted mean across target aggregation areas      
    true_indicators_weighted[,"Head_Count"] <- mapply(FUN=weighted.mean, x=split(pov, pop_domains_vec_tmp),w=split(pop_weights_vec,pop_domains_vec_tmp))
    return(list(true_indicators=true_indicators_weighted,vu_tmp = vu_tmp))
}
# close function 
}


# Superpopulation function -----------------------------------------------------

# The model parameter from the nested error linear regression model are
# used to contruct a superpopulation model.
superpopulation_wild <- function(framework, model_par, gen_model, lambda,
                                 shift, transformation, res_s, fitted_s, fixed){
  # rescaling the errors
  res_s <- sqrt(model_par$sigmae2est) * (res_s - mean(res_s)) / sd(res_s)

  # superpopulation random effect
  vu_tmp <- rnorm(framework$N_dom_pop, 0, sqrt(model_par$sigmau2est))
  vu_pop <- rep(vu_tmp, framework$n_pop)

  # income without individual errors
  Y_pop_b <- gen_model$mu_fixed + vu_pop

  #Commented out because knnx.index in FNN is much faster  
  #indexer <- vapply(Y_pop_b,
  #  function(x) {
  #    which.min(abs(x - fitted_s))
  #  },
  #  FUN.VALUE = integer(1)
  #)

  # index r gets the index of the nearest neighbor of fitted_s for every value of Y_pop_b 
  
  indexer <- FNN::knnx.index(data=fitted_s,query=Y_pop_b,k=1)
  
  
  # superpopulation individual errors
  eps <- res_s[indexer]
  #if (is.null(framework$MSE_cluster)) {
  wu <- sample(c(-1, 1), size = length(eps), replace = TRUE)
  eps <- abs(eps) * wu
  #} 
  #else {
  #}

  #  superpopulation income vector
  Y_pop_b <- Y_pop_b + eps

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

superpopulation <- function(framework, model_par, gen_model, lambda, shift,
                            transformation, fixed) {
  #superpopulation cluster effect 
 
    # superpopulation individual errors
    eps <- vector(length = framework$N_pop)
    eps[framework$obs_dom] <- rnorm(
      sum(framework$obs_dom), 0,
      sqrt(model_par$sigmae2est)
    )
    eps[!framework$obs_dom] <- rnorm(
      sum(!framework$obs_dom), 0,
      sqrt(model_par$sigmae2est +
             model_par$sigmau2est)
    )  
 
  
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


superpopulation_2f <- function(framework, model_par, gen_model, lambda, shift,
                            transformation, fixed) {
  
  # superpopulation area random effect
  vu_tmp <- rnorm(framework$N_dom_pop, 0, sqrt(model_par$sigmau2est))
  vu_pop <- rep(vu_tmp, framework$n_pop)
  
  # superpoulation subarea random effect 
  eta_tmp <- rnorm(framework$N_subdom_pop, 0, sqrt(model_par$sigmah2est))
  eta_pop <- rep(eta_tmp, framework$n_pop_subdom)
  
  # need N_subdom_pop (Number of subdomains in population) and n_pop_subdom (# of units in each subdomain)
  # superpopulation individual errors
  eps <- vector(length = framework$N_pop)
  eps[framework$obs_dom & framework$obs_subdom] <- rnorm(
    sum(framework$obs_dom & framework$obs_subdom), 0,
    sqrt(model_par$sigmae2est)
  )
  eps[framework$obs_dom & !framework$obs_subdom] <- rnorm(
    sum(framework$obs_dom & !framework$obs_subdom), 0,
    sqrt(model_par$sigmae2est+model_par$sigmah2est)
  )
  eps[!framework$obs_dom] <- rnorm(
    sum(!framework$obs_dom), 0,
    sqrt(model_par$sigmae2est +
           model_par$sigmah2est + 
           model_par$sigmau2est)
  )  
  
  
  
  #  superpopulation income vector
  Y_pop_b <- gen_model$mu_fixed + eps + eta_pop + vu_pop
  
  Y_pop_b <- back_transformation(
    y = Y_pop_b,
    transformation = transformation,
    lambda = lambda,
    shift = shift,
    framework = framework,
    fixed = fixed
  )
  Y_pop_b[!is.finite(Y_pop_b)] <- 0
  
  return(list(pop_income_vector = Y_pop_b, vu_tmp = vu_tmp, eta_tmp=eta_tmp))
}



# Bootstrap function -----------------------------------------------------------

bootstrap_par <- function(fixed, transformation, framework, model_par, lambda,
                          shift, vu_tmp) {
  # Bootstrap sample individual error term
  eps <- rnorm(framework$N_smp, 0, sqrt(model_par$sigmae2est))
  # Bootstrap sample random effect
  vu_smp <- rep(vu_tmp[framework$dist_obs_dom], framework$n_smp)
  # Extraction of design matrix
  X_smp <- model.matrix(fixed, framework$smp_data)
  # Constant part of income vector for bootstrap sample
  mu_smp <- X_smp %*% model_par$betas
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




bootstrap_par_wild <- function(fixed, transformation, framework, model_par,
                               lambda, shift, vu_tmp, res_s, fitted_s) {

  # rescaling sample individual error term
  res_s <- sqrt(model_par$sigmae2est) * (res_s - mean(res_s)) / sd(res_s)
  # Bootstrap sample individual error term
  ws <- sample(c(-1, 1), size = length(res_s), replace = TRUE)
  eps <- abs(res_s) * ws

  # Bootstrap sample random effect
  vu_smp <- rep(vu_tmp[framework$dist_obs_dom], framework$n_smp)

  # Extraction of design matrix
  X_smp <- model.matrix(fixed, framework$smp_data)

  # Transformed bootstrap income vector
  Y_smp_b <- X_smp %*% model_par$betas + eps + vu_smp

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

# progress for mse_estim (only internal) ----------

mse_estim_wrapper <- function(i,
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
                              boot_type,
                              seedvec,
                              benchmark,
                              benchmark_type,
                              benchmark_level) {

  tmp <- mse_estim(
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
    boot_type = boot_type,
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
