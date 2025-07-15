# Point estimation function

# This function implements the transformation of data, estimation of the unconditional random effects model 
# and the monte-carlo approximation to predict
# the desired indicators. If the weighted version of the approach is used, then
# additional estimation steps are taken in order to calculate weighted
# regression coefficients before the monte-carlo approximation. See
# corresponding functions below.


point_estim_ell <- function(framework,
                        fixed,
                        alpha, 
                        transformation,
                        interval,
                        L,
                        keep_data = FALSE,
                        Ydump = NULL) {
  
  # Transformation of data -----------------------------------------------------
  
  # Estimating the optimal parameter by optimization
  # Optimal parameter function returns the minimum of the optimization
  # functions from generic_opt; the minimum is the optimal lambda.
  # The function can be found in the script optimal_parameter.R
  optimal_lambda <- optimal_parameter(
    generic_opt = generic_opt,
    fixed = fixed,
    smp_data = framework$smp_data,
    smp_domains = framework$smp_domains,
    transformation = transformation,
    interval = interval,
    framework = framework
  )
  
  # Data_transformation function returns transformed data and shift parameter.
  # The function can be found in the script transformation_functions.R
  transformation_par <- data_transformation(
    fixed = fixed,
    smp_data = framework$smp_data,
    transformation = transformation,
    lambda = optimal_lambda
  )
  shift_par <- transformation_par$shift
   
  
  # Make functions for different estimation_methods 
  
  # Model estimation, model parameter and parameter of generating model --------

  
 
  
if (framework$estimation_method=="plm") {  
  re_model <- estimate_model_plm(formula=fixed,data=transformation_par$transformed_data,
                                 framework=framework)
  est_par <- model_par_ell (
    re_model = re_model,
    framework = framework,
    fixed = fixed,
    transformation_par = transformation_par
  )
  
}
  else if (framework$estimation_method=="gls") {
  gls_model <- estimate_model_gls(formula=fixed,data=transformation_par$transformed_data,
                                  framework=framework)    
  
  est_par <- model_par_ell_gls (
    gls_model = gls_model,
    framework = framework,
    fixed=fixed,
    transformation_par = transformation_par 
  )
  re_model <- gls_model 
  
  }
    

  # Function model_par extracts the needed parameters theta from the random
  # effects linear regression model. It returns the beta coefficients (betas),
  # sigmae2est, and sigmau2est.
  
  
  
  if (!is.null(alpha)) {
    alpha_model <- alphamodel(residuals = est_par$residuals,
                    alpha = alpha,framework = framework)
  }
  else {
    alpha_model <- NULL 
  }
  
  
  
  # Monte-Carlo approximation --------------------------------------------------
  if (inherits(framework$threshold, "function")) {
    framework$threshold <-
      framework$threshold(
        y =
          as.numeric(framework$smp_data[[paste0(fixed[2])]])
      )
  }
  
  # Monte Carlo ELL 
  
  # The monte-carlo function returns a data frame of desired indicators.
  if (L>0) {
    indicator_prediction <- monte_carlo_ell(
      transformation = transformation,
      L = L,
      framework = framework,
      lambda = optimal_lambda,
      shift = shift_par,
      model_par = est_par,
      fixed = fixed,
      Ydump = Ydump,
      alpha_model = alpha_model
    )
  }
  else {
      stop(strwrap(prefix = " ", initial = "",
                   "L must be positive when using ELL estimation"))
    }
  
  
  
  return(list(
    ind = indicator_prediction$point_estimates,
    var = indicator_prediction$var_estimates, 
    optimal_lambda = optimal_lambda,
    shift_par = shift_par,
    model_par = est_par,
    gen_model = gen_model,
    model = re_model,
    alpha_model = alpha_model$alpha_model
    
  ))
} # End point estimation function





# All following functions are only internal ------------------------------------

# functions to estimate random effects model 

estimate_model_plm <- function(formula=formula,data=data,framework=framework) {
  # estimate an unconditional random effect model using plm 
  random_arg <- NULL 
  random_arg[framework$smp_domains] <- list(as.formula(~1))
  names(random_arg) <- c(framework$smp_domains)
  weights_arg <- NULL
  if (!is.null(framework$weight)) {
    weights_arg <- framework$smp_data[,framework$weights]
  }
  
  # Using do.call passes the name of the weight vector from framework to the plm function 
  # if weight is NULL, that is appropriately passed to PLM 
  args <- list(formula=formula, 
               data = data, 
               weights = weights_arg,
               model="random",
               index = framework$smp_domains,
               random.method=framework$random_method)
  
  
  re_model <- do.call(plm:::plm, args)
  
  
  return(re_model)
} # end function estimate_model_plm 


estimate_model_gls <- function(formula=formula,data=data,framework=framework) {
  # estimate an unconditional random effects model using gls 
  weights_arg <- NULL
  if (!is.null(framework$weight)) {
    data$weights_temp <- framework$smp_data[,framework$weights]
    weights_arg <- ~1/weights_temp
  }
  correlation_form <- as.formula(paste0("~1 |",framework$smp_domains))
  correlation_arg <- paste0("nlme:::corCompSymm(form= ~ 1 | TA_CODE)")
  coefs <- all.vars(formula)
  summary(data[,coefs])
  
  ## lme_gls <- nlme:::gls(model=ebpformula,data=estimationSample,correlation =  corCompSymm(form = ~ 1 | TA_CODE),method="ML")
  
  corr_form <- as.formula(paste0("~1 | ",framework$smp_domains))
  corr_arg <- nlme:::corCompSymm(form= corr_form)
  
  args <- list(model=formula,
               data=data,
               weights = weights_arg, 
               method = "ML",
               correlation = corr_arg
  ) 
  gls_model <- do.call(nlme:::gls,args)
  return(gls_model)
} # end function estimate_model_gls 




# Functions to extract and calculate model parameter----------------------------

# Function model_par extracts the needed parameters theta from the nested
# error linear regression model. It returns the beta coefficients (betas),
# sigmae2est, sigmau2est and the random effect (rand_eff).

model_par_ell <- function(framework,
                      re_model,
                      fixed,
                      transformation_par) {
  
  # fixed parameters
  betas <- re_model$coefficients 
  # Estimated error variance
  sigmae2est <- re_model$ercomp$sigma2[1]
  # VarCorr(fit2) is the estimated random error variance
  sigmau2est <- re_model$ercomp$sigma2[2]
  
  if (framework$model_parameters!="fixed") {
    varFix <- re_model$vcov
  }
  else {
    varFix = NULL 
  }
  residuals <- as.vector(re_model$model[,1]-predict(re_model))
  if (is.null(framework$weights)) {
    weights <- rep(1,framework$N_smp)
  }
  else {
    weights <- framework$smp_data[,framework$weights]
  }
  mean_residuals <- aggregate_weighted_mean(df=residuals,by=list(framework$smp_data[,framework$smp_domains]),
                                        w=weights)[,2]
  framework$pop_data[[paste0(fixed[2])]] <- seq_len(nrow(framework$pop_data))
  X_pop <- model.matrix(fixed, framework$pop_data)
  mu_fixed <- X_pop %*% betas
  
  
    return(list(
      betas = betas,
      sigmae2est = sigmae2est,
      sigmau2est = sigmau2est,
      varFix = varFix,
      varErr = NULL,
      residuals = residuals,
      mean_residuals = mean_residuals, 
      mu_fixed = mu_fixed 
    ))
  } 

model_par_ell_gls <- function(framework,
                          gls_model,
                          fixed,
                          transformation_par) {
  
  # fixed parameters
  betas <- gls_model$coefficients
  # variance parameters 
  rho<-coef(lme_gls$modelStruct,unconstrained=F)
  # area effect 
  sigmau2est <- lme_gls$sigma^2*rho
  # Estimated epsilon error variance
  sigmae2est <- lme_gls$sigma^2*(1-rho)
  
  if (framework$model_parameters!="fixed") {
    varFix <- gls_model$varBeta
  }
  else {
    varFix = NULL 
  }
  residuals <- gls_model$residuals
  if (is.null(framework$weights)) {
    weights <- rep(1,framework$N_smp)
  }
  else {
    weights <- framework$smp_data[,framework$weights]
  }
  mean_residuals <- aggregate_weighted_mean(df=residuals,by=list(framework$smp_data[,framework$smp_domains]),
                                            w=weights)[,2]
  framework$pop_data[[paste0(fixed[2])]] <- seq_len(nrow(framework$pop_data))
  X_pop <- model.matrix(fixed, framework$pop_data)
  mu_fixed <- X_pop %*% betas
  
  
  return(list(
    betas = betas,
    sigmae2est = sigmae2est,
    sigmau2est = sigmau2est,
    varFix = varFix,
    varErr = NULL,
    residuals = residuals,
    mean_residuals = mean_residuals, 
    mu_fixed = mu_fixed 
  ))
} 






# function rowvar 
# calcluate the variance of a row 
rowvar <- function(x) {
  rowvar <-apply(x,1,var)
}


#Alpha model function
# This function estimates the alpha model, as described in Zhao (2006) p.4-5 
#Zhao, Q. (2006). User manual for povmap. World Bank. http://siteresources. worldbank. org/INTPGI/Resources/342674-1092157888460/Zhao_ ManualPovMap. pdf.
# it returns the alpha model and the expected value of the variance 
alphamodel <- function(residuals, alpha,framework) {
  if (is.null(framework$weights)) {
    weight_smp <- rep(1,framework$N_smp)  
  } 
  else {
    weight_smp <- framework$smp_data[[as.character(framework$weights)]]
  }
  
    # 1. Decompose the residuals into an average cluster effect and a residual 
  
   mean_resid <- aggregate_weighted_mean(df=residuals,by=list(framework$smp_data[,framework$smp_domains]),
                                         w=weight_smp)
   dev_resid <- residuals - rep(mean_resid$V1,framework$n_smp)
   
   eps_squared <- dev_resid^2 
    A <- 1.05*max(eps_squared)
    framework$smp_data$transformed_eps_squared <- log(eps_squared/(A-eps_squared))
    framework$smp_data$transformed_eps_squared[eps_squared==0] <- 0
    model <- reformulate(response="transformed_eps_squared",as.character(alpha[-1]))
    
    alpha_model<-lm(model,data=framework$smp_data,weights=weight_smp)
    
    # we want to draw standardized residuals, so first estimate the variance of epsilon in the sample 
    #B is not explicitly defined in Zhang (2006) but Peter confirmed it is exp(zTch*alphahat) as defined in Zhang (2006)
    # and that has to be true for the variance formula to be correct 
    B_smp <- exp(predict(alpha_model))
    
    var_r <- summary(alpha_model)$sigma^2
    # This formula comes from Zhao (2006): User Manual for Povmap 
    sigmae2est_smp <- (A * B_smp / (1+B_smp)) + 0.5*var_r*(A*B_smp*(1-B_smp)/(1+B_smp)^3)
    dev_resid_std <- dev_resid/sigmae2est_smp^0.5
    dev_resid_std <- dev_resid_std-ave(dev_resid_std,by=framework$smp_data[,framework$smp_domains]) # ELL p.357 eq (3)
      
    # now we want to unstandardize the residuals, so we need the estimated variance of epsilon in the population
    alpha_X_vars <- alpha_model$terms
    framework$pop_data[[paste0(alpha_X_vars[[2]])]] <- seq_len(nrow(framework$pop_data))
    X_pop <- model.matrix(alpha_X_vars, framework$pop_data)
    B_pop <- exp(X_pop %*% alpha_model$coefficients)
    sigmae2est_pop <- (A * B_pop / (1+B_pop)) + 0.5*var_r*(A*B_pop*(1-B_pop)/(1+B_pop)^3)
    sigmae2est_pop[sigmae2est_pop<min(sigmae2est_smp)] <- min(sigmae2est_smp)
    sigmae2est_pop[sigmae2est_pop>max(sigmae2est_smp)] <- max(sigmae2est_smp)
    return(list(alpha_model=alpha_model,sigmae2est_smp=sigmae2est_smp, sigmae2est_pop=sigmae2est_pop, 
                dev_resid_std=dev_resid_std,mean_resid = mean_resid$V1))
    }


# Monte-Carlo approximation ----------------------------------------------------

# This function conducts bootstraps 
monte_carlo_ell <- function(transformation,
                        L,
                        framework,
                        lambda,
                        shift,
                        model_par,
                        fixed,
                        Ydump,
                        alpha_model = NULL) {
  
  # Preparing matrices for indicators for the Monte-Carlo simulation
  
  if(!is.null(framework$aggregate_to_vec)){
    N_dom_pop_tmp <- framework$N_dom_pop_agg
    pop_domains_vec_tmp <- framework$aggregate_to_vec
  } else {
    N_dom_pop_tmp <- framework$N_dom_pop
    pop_domains_vec_tmp <- framework$pop_domains_vec
  }
  
  
  if (!is.null(Ydump)) {
    Ydumpdf <- data.frame(matrix(ncol = 6, nrow = 0))
    colnames(Ydumpdf) <- c("L","Domain","Simulated_Y","XBetahat","eta","epsilon")
    write.csv(Ydumpdf,Ydump,row.names = FALSE)
  }
  
  
  ests_mcmc <- array(dim = c(
    N_dom_pop_tmp,
    L,
    length(framework$indicator_names)
  ))
  
  betas <- model_par$betas 
  message("\r", "Bootstrap started                                            ")
  start_time <- Sys.time()
  for (l in seq_len(L)) {
    
    if (framework$model_parameters=="variable") {
      # variable parameters means they must be drawn every replication
      # so we redraw the parameters  
      # draw error terms 
      #R <- chol(model_par$varErr) 
      #sigma2<- c(-1,-1)
      #while (sigma2[2]<0 | sigma2[1]<0) {
      #  sigma2 <- t(R)  %*% matrix(rnorm(ncol(R)), ncol(R)) + diag(model_par$varErr)
      #  model_par$sigmae2est <- sigma2[2]
      #  model_par$sigmau2est <- sigma2[1]
      #}
      # draw betas 
      R <- chol(model_par$varFix)
      model_par$betas <- t(R)  %*% matrix(rnorm(ncol(R)), ncol(R)) + betas  
      
    } # close condition to redraw parameters 
    
    
    
    # Errors in generating model: individual error term and random effect
    # See below for function errors_gen_ell.
    
    if (framework$errors!="nonnormal") {
    errors <- errors_gen_ell(
      framework = framework,
      model_par = model_par,
      alpha_model = alpha_model
    )
    } 
    else {
      errors <- errors_gen_ell_nonp(
        framework=framework,
        model_par = model_par,
        alpha_model = alpha_model 
      )
    
    }
    
    # generate gen_model$mu here based parameters$beta, if you have the data conveniently around 
    
    
    
    # Prediction of population vector y
    # See below for function prediction_y_ell.
    population_vector <- prediction_y_ell(
      transformation = transformation,
      lambda = lambda,
      shift = shift,
      model_par = model_par, 
      errors = errors,
      framework = framework,
      fixed = fixed
    )
    
    if(!is.null(framework$pop_weights)){
      pop_weights_vec <- framework$pop_data[[framework$pop_weights]]
    }else{
      pop_weights_vec <- rep(1, nrow(framework$pop_data))
    }
    if (!is.null(Ydump)){
      Ydumpdf <- data.frame(rep(l,framework$N_pop), framework$pop_domains_vec,population_vector,model_par$mu_fixed,errors$vu,errors$epsilon)
      #write.csv(Ydumpdf,Ydump,row.names = FALSE,append=TRUE)
      write.table(Ydumpdf,file=Ydump,row.names = FALSE,append=TRUE,col.names=F, sep=",") 
    }
    
    # Calculation of indicators for each Monte Carlo population
    ests_mcmc[, l, ] <-
      matrix(
        nrow = N_dom_pop_tmp,
        data = unlist(lapply(framework$indicator_list,
                             function(f, threshold) {
                               matrix(
                                 nrow = N_dom_pop_tmp,
                                 data = unlist(mapply(
                                   y = split(population_vector, pop_domains_vec_tmp),
                                   pop_weights = split(pop_weights_vec, pop_domains_vec_tmp),
                                   f,
                                   threshold = framework$threshold
                                 )), byrow = TRUE
                               )
                             },
                             threshold = framework$threshold
        ))
      )
  
    
    if (l %% 20 == 0) {
      if (l != L) {
        delta <- difftime(Sys.time(), start_time, units = "secs")
        remaining <- (delta / l) * (L - l)
        remaining <- unclass(remaining)
        remaining <- sprintf(
          "%02d:%02d:%02d:%02d",
          remaining %/% 86400, # days
          remaining %% 86400 %/% 3600, # hours
          remaining %% 3600 %/% 60, # minutes
          remaining %% 60 %/% 1
        ) # seconds)
        
        message("\r", l, " of ", L, " Bootstrap iterations completed \t
              Approximately ", remaining, " remaining \n")
        #if (.Platform$OS.type == "windows") flush.console()
      }
    }
    
  } # End for loop
  
  # Point estimations of indicators by taking the mean
  
  point_estimates <- data.frame(
    Domain = unique(pop_domains_vec_tmp),
    apply(ests_mcmc, c(3), rowMeans)
  )
  var_estimates <- data.frame(
    Domain=unique(pop_domains_vec_tmp),
    apply(ests_mcmc,3,rowvar))
  
  
  colnames(point_estimates) <- c("Domain", framework$indicator_names)
  colnames(var_estimates) <- c("Domain", framework$indicator_names)
  

  
  return(list(point_estimates = point_estimates, var_estimates = var_estimates))
} # End Monte-Carlo


# The function errors_gen returns error terms of the generating model.
errors_gen_ell <- function(framework, model_par, alpha_model=NULL) {
  if (is.null(alpha_model)) {
    epsilon <- rnorm(framework$N_pop, 0, sqrt(model_par$sigmae2est))
  } 
  else {
    epsilon <- rnorm(framework$N_pop, 0, sqrt(alpha_model$sigmae2est_pop))
  }
  # empty vector for new random effect in generating model
  vu <- vector(length = framework$N_pop)
  # draw random effect
  vu <- rep(rnorm(framework$N_dom_pop,0,sqrt(model_par$sigmau2est)
    ),framework$n_pop)
  return(list(epsilon = epsilon, vu = vu))
} # End errors_gen_ell

# The function errors_gen_nonp returns error terms of the generating model
# obtained through a non-parametric bootstrap procedure 
errors_gen_ell_nonp <- function(framework, model_par, alpha_model) {
 
  if (!is.null(alpha_model)) {
    epsilon_std <- sample(alpha_model$dev_resid_std, replace=TRUE,size=framework$N_pop)
    # now we want to unstandardize the residuals, so we need the estimated variance of epsilon in the population
    epsilon <- as.vector(epsilon_std*alpha_model$sigmae2est_pop^0.5)
    vu <- rep(sample(alpha_model$mean_resid,replace=TRUE,size=framework$N_dom_pop),framework$n_pop) 
  }
     else {
       weights <- framework$smp_data[,framework$weights]
       residuals <- model_par$residuals 
       mean_resid <- aggregate_weighted_mean(df=residuals,by=list(framework$smp_data[,framework$smp_domains]),
                                             w=weights)
       dev_resid <- residuals - rep(mean_resid$V1,framework$n_smp)
       epsilon <- sample(dev_resid, replace=TRUE,size=framework$N_pop)
       vu <- rep(sample(mean_resid$V1,replace=TRUE,size=framework$N_dom_pop),framework$n_pop) 
     }
   
  return(list(epsilon = epsilon, vu = vu))
} # End errors_gen_ell_nonp





# The function prediction_y returns a predicted income vector which can be used
# to calculate indicators. Note that a whole income vector is predicted without
# distinction between in- and out-of-sample domains.
prediction_y_ell <- function(transformation,
                         lambda,
                         shift,
                         model_par,
                         errors,
                         framework,
                         fixed) {
  # construct vector of Xes by starting with intercept and binding it to the population
  # dataframe constructed from all the coefficient names except the 

  
  # predicted population income vector
  y_pred <- model_par$mu_fixed + errors$epsilon + errors$vu
  
  # back-transformation of predicted population income vector
  y_pred <- back_transformation(
    y = y_pred,
    transformation = transformation,
    lambda = lambda,
    shift = shift,
    framework = framework,
    fixed = fixed
  )
  y_pred[!is.finite(y_pred)] <- 0
  
  return(y_pred)
} # End prediction_y_ell
