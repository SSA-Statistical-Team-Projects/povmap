# Internal documentation -------------------------------------------------------

#quiet <- function(x) { 
#  sink(tempfile()) 
#  on.exit(sink()) 
#  invisible(force(x)) 
#} 

# Point estimation function

# This function implements the transformation of data, estimation of the nested
# error linear regression model and the monte-carlo approximation to predict
# the desired indicators. If the weighted version of the approach is used, then
# additional estimation steps are taken in order to calculate weighted
# regression coefficients before the monte-carlo approximation. See
# corresponding functions below.


point_estim <- function(framework,
                        fixed,
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

  # Model estimation, model parameter and parameter of generating model --------

  # Estimation of the nested error linear regression model
  # See Molina and Rao (2010) p. 374
  # lme function is included in the nlme package which is imported.
  weights_arg <- NULL 
  if(!is.null(framework$weights) &&
     any(framework$weights_type %in% c("nlme", "nlme_lambda","hybrid2"))) {
    #transformation_par$transformed_data$weights_scaled <-
    #  framework$smp_data[,framework$weights] /
    #    mean(framework$smp_data[,framework$weights], na.rm = TRUE)
    transformation_par$transformed_data$weights_temp <- framework$smp_data[,framework$weights]
    #weights_arg <- nlme:::varComb(nlme:::varIdent(form = ~1),nlme:::varFixed(~1/weights_temp))
    weights_arg <- ~1/weights_temp
  }   
  
  random_arg <- NULL 
  if (!is.null(framework$smp_subdomains) && !is.null(framework$pop_subdomains)) {
    # Do two fold model 
    random_arg <- list(as.formula(~1),as.formula(~1))
    names(random_arg) <- c(framework$smp_domains,framework$smp_subdomains)  
    
      }
  else {
  # Do one-fold model 
  random_arg[framework$smp_domains] <- list(as.formula(~1))
  } 
  
  
  args <- list(fixed=fixed,
               data = transformation_par$transformed_data,
               random = random_arg,
               method = framework$nlme_method,
               control = nlme::lmeControl(maxIter = framework$nlme_maxiter,
                                          tolerance = framework$nlme_tolerance,
                                          opt = framework$nlme_opt,
                                          optimMethod = framework$nlme_optimmethod, 
                                          msMaxIter=framework$nlme_msmaxiter,
                                          msTol=framework$nlme_mstol,
                                          returnObject = framework$nlme_returnobject 
               ),
               keep.data = keep_data,
               weights = weights_arg)
  
  #ptm <- proc.time()
  mixed_model <- do.call(nlme:::lme,args)
  #proc.time() - ptm     
 
    
    # Function model_par extracts the needed parameters theta from the nested
    # error linear regression model. It returns the beta coefficients (betas),
    # sigmae2est, sigmau2est, sigmah2est (which is set to 0 for one-fold models) and the random effect (rand_eff)

    dep_var <- transformation_par$transformed_data[[
      as.character(mixed_model$terms[[2]])]]
       # this is needed for guadarrama weights calculation in gen_model 


  # Function gen_model calculates the parameters in the generating model.
  # See Molina and Rao (2010) p. 375 (20)
  # The function returns sigmav2est and the constant part mu.
    est_par <- model_par(
      mixed_model = mixed_model,
      framework = framework,
      fixed = fixed,
      transformation_par = transformation_par
    )
    
  if (framework$model_parameters!="variable") {
    gen_par <- gen_model(
    model_par = est_par,
    fixed = fixed,
    framework = framework,
    dep_var = dep_var
  )
  } 
  else {
    gen_par <- NULL 
  }
  #update random effects, beta coefficients, variance components, and residuals in the cases that they are altered by gen_model
  est_par$rand_eff <- gen_par$rand_eff 
  est_par$betas <- gen_par$betas
  est_par$sigmau2est <- gen_par$sigmau2est
  est_par$sigmae2est <- gen_par$sigmae2est
  est_par$e0 <- gen_par$e0
  est_par$e1 <- gen_par$e1
  est_par$e2 <- gen_par$e2
  est_par$gamma <- gen_par$gamma 
  est_par$gamma_sub <- gen_par$gamma_sub
  
  
  # Monte-Carlo approximation --------------------------------------------------
  if (inherits(framework$threshold, "function")) {
    framework$threshold <-
      framework$threshold(
        y =
          as.numeric(framework$smp_data[[paste0(fixed[2])]])
      )
  }

  # The monte-carlo function returns a data frame of desired indicators.
 if (L>0 & framework$data.table==FALSE) {
   indicator_prediction <- monte_carlo(
    mixed_model = mixed_model, 
    transformation = transformation,
    transformation_par = transformation_par, 
    L = L,
    framework = framework,
    lambda = optimal_lambda,
    shift = shift_par,
    model_par = est_par,
    gen_model = gen_par,
    fixed = fixed,
    Ydump = Ydump 
  )
 }
  else if (L>0 & framework$data.table==TRUE) {
    indicator_prediction <- monte_carlo_dt(
      mixed_model=mixed_model,
      transformation = transformation,
      transformation_par = transformation_par, 
      L = L,
      framework = framework,
      lambda = optimal_lambda,
      shift = shift_par,
      model_par = est_par,
      gen_model = gen_par,
      fixed = fixed,
      Ydump = Ydump 
    )
    
    
  }
    else {
    indicator_prediction <- analytic(
      transformation=transformation,
      framework = framework,
      lambda = optimal_lambda,
      shift = shift_par,
      model_par = est_par,
      gen_model = gen_par, 
      fixed = fixed, 
      Ydump=Ydump
    )
  }
  
  
  mixed_model$coefficients_weighted <- if (!is.null(framework$weights)) {
    as.numeric(est_par$betas)
  } else {
    NULL
  }
  names(mixed_model$coefficients_weighted) <- if (!is.null(framework$weights)) {
    rownames(est_par$betas)
  } else {
    NULL
  }
 
  
  return(list(
    ind = indicator_prediction,
    optimal_lambda = optimal_lambda,
    shift_par = shift_par,
    model_par = est_par,
    gen_model = gen_par,
    model = mixed_model
  ))
} # End point estimation function


# All following functions are only internal ------------------------------------

# Functions to extract and calculate model parameter----------------------------

# Function model_par extracts the needed parameters theta from the nested
# error linear regression model. It returns the beta coefficients (betas),
# sigmae2est, sigmau2est and the random effect (rand_eff).

model_par <- function(framework,
                      mixed_model,
                      fixed,
                      transformation_par) {

  # fixed parameters
  betas <- nlme::fixed.effects(mixed_model)
  # Estimated error variance
  sigmae2est <- mixed_model$sigma^2
  # VarCorr(fit2) is the estimated random error variance
  # if subdomains specified, save the results in sigmae2est, sigmah2est, and sigmau2est 
  # Random effect: vector with zeros for all domains, filled with 0
  rand_eff <- rep(0, framework$N_dom_pop)  
  rand_eff_h <- rep(0, framework$N_subdom_pop)  
  rand_eff_h_smp <- NULL 
  if (!is.null(framework$smp_subdomains) && !is.null(framework$pop_subdomains)) {
    sigmau2est <- as.numeric(VarCorr(mixed_model)[2,1])
    sigmah2est <- as.numeric(VarCorr(mixed_model)[4,1])
    # Random area effect for sample domains
    rand_eff[framework$dist_obs_dom] <- random.effects(mixed_model,level=1)$"(Intercept)"[framework$dist_obs_smp_dom]
    # Random subarea-effect 
    rand_eff_h[framework$dist_obs_subdom] <- random.effects(mixed_model,level=2)$"(Intercept)"[framework$dist_obs_smp_subdom]
    rand_eff_h_smp <- random.effects(mixed_model,level=2)$"(Intercept)"
  } 
  else {
    sigmau2est <- as.numeric(nlme::VarCorr(mixed_model)[1, 1])  
    sigmah2est <- 0 
    # random area effect 
    rand_eff[framework$dist_obs_dom] <- (random.effects(mixed_model)[[1]])
  }
  if (framework$model_parameters!="fixed") {
    varFix=mixed_model$varFix 
    varErr = lmeInfo:::varcomp_vcov(mixed_model)
      }
  else {
    varFix = NULL 
    varErr = NULL 
  }
  
 e0 <- residuals(mixed_model,level=0) 
 
 
  
  
    
    
    return(list(
      betas = betas,
      sigmae2est = sigmae2est,
      sigmah2est = sigmah2est, 
      sigmau2est = sigmau2est,
      rand_eff = rand_eff,
      rand_eff_h = rand_eff_h,
      rand_eff_h_smp = rand_eff_h_smp, 
      varFix=varFix, 
      varErr=varErr,
      e0=e0
    ))

} # End model_par



  
  
  
# Function gen_model calculates the parameters in the generating model with a one-fold model.
# See Molina and Rao (2010) p. 375 (20)
gen_model <- function(fixed,
                      framework,
                      model_par,
                      dep_var) {

  if (is.null(framework$weights)) {
    weight_smp <- rep(1,framework$N_smp)  
  } 
  else {
    weight_smp <- framework$smp_data[[as.character(framework$weights)]]
  }
    
  betas <- model_par$betas
  sigmae2est <- model_par$sigmae2est 
  sigmau2est <- model_par$sigmau2est 
  rand_eff <- model_par$rand_eff
  rand_eff_h <- model_par$rand_eff_h 
  gamma_sub <- NULL 
  
  # calculate gamma 
  if (any(framework$weights_type %in% c("nlme", "nlme_lambda","hybrid2",NULL)) | is.null(framework$weights)) {
  weight_sum <- rep(0, framework$N_dom_smp)
  sums <- aggregate(data.frame(weight_smp,weight_smp^2), by=list(framework$smp_domains_vec),FUN=sum)
  delta2 <- sums[,3] / sums[,2]^2 # sum of the squares divided by the square of the sum
  #when using wrong_gamma, random effects do not change, showing that nlme by default ignores weights when calculating gamma 
  #ones<-rep(1,length(weight_smp))
  #wrong_sums <- aggregate(data.frame(ones,ones), by=list(framework$smp_domains_vec),FUN=sum)
  #wrong_delta2 <- wrong_sums[,3] / wrong_sums[,2]^2
  #wrong_gamma <- model_par$sigmau2est / (model_par$sigmau2est + ((model_par$sigmae2est + model_par$sigmah2est) * wrong_delta2))
  gamma <- model_par$sigmau2est / (model_par$sigmau2est + ((model_par$sigmae2est + model_par$sigmah2est) * delta2))
  
  #gamma <- model_par$sigmau2est / (model_par$sigmau2est + ((model_par$sigmae2est + model_par$sigmah2est) * wrong_delta2))
  if (model_par$sigmah2est>0) {
    sums_sub <- aggregate(data.frame(weight_smp, weight_smp^2)[framework$smp_subdomains_vec,], by = list(framework$smp_subdomains_vec), FUN = sum)
    sums_sub <- sums_sub[framework$dist_obs_smp_subdom,]
    delta2_sub <- sums_sub[,3] / sums_sub[,2]^2
    gamma_sub <- model_par$sigmah2est / (model_par$sigmah2est + model_par$sigmae2est * delta2_sub)
  }
   if (!is.null(framework$weights_type)) {
    if (framework$weights_type=="hybrid2") {
      #First update betas 
      newbetas <- update_beta_re(model_par = model_par, weight_smp=weight_smp, framework=framework,dep_var=dep_var,fixed=fixed)
      betas <- newbetas$betas
      gamma <- newbetas$gamma
      rand_eff <- newbetas$rand_eff
      rand_eff_h <- newbetas$rand_eff_h 
      # now update variance components 
      updated_sigma2 <- update_sigma2(dep_var=dep_var,betas=betas,weight_smp=weight_smp,framework=framework,rand_eff=rand_eff,fixed=fixed)
      sigmau2est <- updated_sigma2$sigmau2est 
      sigmae2est <- updated_sigma2$sigmae2est 
    } # close hybrid2 
   } # close non-null 
} # close nlme family 
  else if (!is.null(framework$weights_type)) {
  if (framework$weights_type=="hybrid") {
    newbetas <- update_beta_re(model_par = model_par, weight_smp=weight_smp, framework=framework,dep_var=dep_var,fixed=fixed)
    betas <- newbetas$betas
    gamma <- newbetas$gamma 
    rand_eff <- newbetas$rand_eff
    rand_eff_h <- newbetas$rand_eff_h 
    updated_sigma2 <- update_sigma2(dep_var=dep_var,betas=betas,weight_smp=weight_smp,framework=framework,rand_eff=rand_eff,fixed=fixed)
    sigmau2est <- updated_sigma2$sigmau2est 
    sigmae2est <- updated_sigma2$sigmae2est 
    } # close hybrid
    else if (framework$weights_type=="Guadarrama") {
      # Calculations needed for pseudo EB for Guadarrama option 
      weight_sum <- rep(0, framework$N_dom_smp)
      mean_dep <- rep(0, framework$N_dom_smp)
      mean_indep <- matrix(0, nrow = framework$N_dom_smp, ncol = length(betas))
      delta2 <- rep(0, framework$N_dom_smp)
      gamma <- rep(0, framework$N_dom_smp)
      num <- matrix(0, nrow = length(betas), ncol = 1) 
      den <- matrix(0, nrow = length(betas), ncol = length(betas))
      
      for (d in 1:framework$N_dom_smp) {
        #cat(d)
        domain <- names(table(framework$smp_domains_vec)[d])

        # Domain means of of the dependent variable
        dep_smp <- dep_var[framework$smp_domains_vec == domain]
        weight_smp <- framework$smp_data[[
          as.character(framework$weights)]][framework$smp_domains_vec == domain]
        weight_sum[d] <- sum(weight_smp)
        
         if(length(weight_smp) == 1) {
           indep_smp <- matrix(model.matrix(fixed, framework$smp_data)[framework$smp_domains_vec == domain,]
                               , ncol = length(betas), nrow = 1)
        } else {
          indep_smp <- model.matrix(fixed, framework$smp_data)[framework$smp_domains_vec == domain,]
        }
        
        # weighted mean of the dependent variable
        mean_dep[d] <- sum(weight_smp * dep_smp) / weight_sum[d]
        
        # weighted means of the auxiliary information
        for (k in 1:length(betas)) {
          mean_indep[d, k] <- sum(weight_smp * indep_smp[, k]) / weight_sum[d]
        }
        
        delta2[d] <- sum(weight_smp^2) / (weight_sum[d]^2)
        gamma[d] <- model_par$sigmau2est / (model_par$sigmau2est + model_par$sigmae2est * delta2[d])
        weight_smp_diag <- diag(x=weight_smp,nrow=length(weight_smp), ncol=length(weight_smp))
        dep_var_ast <- dep_smp - gamma[d] * mean_dep[d]
        indep_weight <- t(indep_smp) %*% weight_smp_diag
        
        a <- matrix(rep(
          gamma[d] *
            mean_indep[d, ],
          framework$n_smp[d]
        ),
        nrow = framework$n_smp[d],
        byrow = TRUE
        )
        
        indep_var_ast <- indep_smp - matrix(rep(
          gamma[d] *
            mean_indep[d, ],
          framework$n_smp[d]
        ),
        nrow = framework$n_smp[d],
        byrow = TRUE
        )
        
        
        
        num <- num + (indep_weight %*% dep_var_ast)
        den <- den + (indep_weight %*% indep_var_ast)
      } # close loop over domains 
      
      
      betas <- solve(den) %*% num
      
      # random effect for in-sample domains (dist_obs_dom)
      rand_eff <- rep(0, framework$N_dom_pop)  
      rand_eff[framework$dist_obs_dom] <- gamma * (mean_dep -
                                                            mean_indep %*% betas)
      
    
    
  } # close Guadarrama   
  } # close non-null

  framework$pop_data[[paste0(fixed[2])]] <- seq_len(nrow(framework$pop_data))
  X_pop <- model.matrix(fixed, framework$pop_data)
  
  # Constant part of predicted y
  mu_fixed <- X_pop %*% betas
  sigmav2est <- sigmau2est * (1 - gamma)
  rand_eff_pop <- rep(rand_eff, framework$n_pop)
  X_smp <- model.matrix(fixed, framework$smp_data)
  e0 <- dep_var - X_smp %*% betas 
  e2 <- NULL  
  
  if (model_par$sigmah2est==0) {
    mu <- mu_fixed + rand_eff_pop
    sigmai2est <- 0 
    e1 <- e0 - rep(rand_eff[framework$dist_obs_dom],framework$n_smp)
      }
  else {
    sigmai2est <- model_par$sigmah2est * (1-gamma_sub)
    rand_eff_h_pop <- rep(rand_eff_h,framework$n_pop_subdom)
    mu <- mu_fixed + rand_eff_pop + rand_eff_h_pop 
    e1 <- e0 - rep(rand_eff[framework$dist_obs_dom],framework$n_smp)
    #e2 <- e1 - rep(rand_eff_h[framework$dist_obs_subdom],framework$n_smp_subdom)
    e2 <- e1 - rep(model_par$rand_eff_h_smp,framework$n_smp_subdom)
      }

  
 
  
  
      return(list(betas=betas,sigmau2est=sigmau2est,sigmae2est=sigmae2est,sigmav2est = sigmav2est, sigmai2est = sigmai2est, 
                  mu = mu, mu_fixed = mu_fixed,rand_eff=rand_eff,gamma=gamma,gamma_sub=gamma_sub,e0=e0,e1=e1,e2=e2))
} # End gen_model


#calculate weighted mean using R's aggregate function 
aggregate_weighted_mean  <-function(df,by,w) {
  wdf <- cbind(df*w,w)  
  wdf_sum <- aggregate(wdf,by=by, FUN=sum)
  aggregate_weighted_mean <- cbind(Domain = wdf_sum[,1],wdf_sum[,-1]/wdf_sum$w)
  aggregate_weighted_mean <-within(aggregate_weighted_mean,rm(w))
  return(aggregate_weighted_mean)
}

# calculated weighted quantile 
aggregate_weighted_quantile  <-function(df,by,w,q) {
   aggregate_weighted_quantile <- mapply(FUN=wtd.quantile,x=split(df,by),weight=split(w,pop_domains_vec_tmp),probs=q)
}
  
#Function to update betas and random effects using Guadarrama et al method 
update_beta_re <- function(model_par,weight_smp,framework,dep_var,fixed) {
  #This code implements Guadarrama et al, starting with parameters from weighted nlme estimated above 
  rand_eff <- model_par$rand_eff
  rand_eff_h <- model_par$rand_eff_h
  sums <- aggregate(data.frame(weight_smp, weight_smp^2), 
                    by = list(framework$smp_domains_vec), FUN = sum)
  delta2 <- sums[, 3]/sums[, 2]^2
  gamma <- model_par$sigmau2est/(model_par$sigmau2est + 
                                   ((model_par$sigmae2est + model_par$sigmah2est) * 
                                      delta2))
  indep_smp <- model.matrix(fixed, framework$smp_data)
  mean_indep <- aggregate_weighted_mean(indep_smp,by=list(framework$smp_domains_vec),w=weight_smp)[,-1]
  # weighted mean of the dependent variable
  mean_dep <- aggregate_weighted_mean(dep_var,by=list(framework$smp_domains_vec),w=weight_smp)[,-1]
  dep_var_ast <- dep_var -  rep(gamma * mean_dep,framework$n_smp)
  #weighted independent variables     
  indep_weight <- indep_smp * weight_smp
  # shrink weighted independent variables 
  shrunk_mean_indep <- gamma*mean_indep 
  # expand from one observation per domain to one observation per sample household 
  shrunk_mean_indep_smp <- shrunk_mean_indep[rep(row.names(shrunk_mean_indep), times = framework$n_smp), ]
  indep_var_ast <- as.matrix(indep_smp - shrunk_mean_indep_smp)  
  
  # If two fold model, subtract subarea gamma from indep_var_ast and dep_va
  if (model_par$sigmah2est>0) {
    sums_sub <- aggregate(data.frame(weight_smp, weight_smp^2)[framework$smp_subdomains_vec,], by = list(framework$smp_subdomains_vec), FUN = sum)
    sums_sub <- sums_sub[framework$dist_obs_smp_subdom,]
    delta2_sub <- sums_sub[,3] / sums_sub[,2]^2
    gamma_sub <- model_par$sigmah2est / (model_par$sigmah2est + model_par$sigmae2est * delta2_sub)
    submean_dep <- aggregate_weighted_mean(dep_var,by=list(framework$smp_subdomains_vec),w=weight_smp)[,-1]
    dep_var_ast <- dep_var_ast - rep(gamma_sub * submean_dep,framework$n_smp_subdom)
    mean_indep_sub <- aggregate_weighted_mean(indep_smp,by=list(framework$smp_subdomains_vec),w=weight_smp)[,-1]
    shrunk_mean_indep_sub <- gamma_sub*mean_indep_sub 
    shrunk_mean_indep_sub_smp <- shrunk_mean_indep_sub[rep(row.names(shrunk_mean_indep_sub), times = framework$n_smp_subdom), ]
    indep_var_ast <- as.matrix(indep_var_ast-shrunk_mean_indep_sub_smp)
  }
  
  num <- t(indep_weight) %*% dep_var_ast
  den <- t(indep_weight) %*% indep_var_ast
  betas <- solve(den) %*% num
  rand_eff[framework$dist_obs_dom] <- gamma * (mean_dep -
                                                 as.matrix(mean_indep) %*% betas)
   # also update random effects for sub-area models 
   
  if (model_par$sigmah2est>0) {
          mean_e0_sub <- aggregate_weighted_mean(model_par$e0,by=list(framework$smp_subdomains_vec),w=weight_smp)
          rand_eff_h[framework$dist_obs_subdom] <- gamma_sub*mean_e0_sub[,2]      
        }
                                               
  return(list(betas=betas,rand_eff=rand_eff,rand_eff_h=rand_eff_h,gamma=gamma))
}

#function to update sigma2 
update_sigma2 <- function(dep_var, betas, weight_smp, framework,rand_eff,fixed) {
indep_smp <- model.matrix(fixed, framework$smp_data)
e0 <- dep_var - indep_smp %*% betas
rand_eff_smp <- rep(rand_eff[framework$dist_obs_dom],framework$n_smp)
#e1 <- e0 - rand_eff_smp  
#e2 = e0 - e1
#w2=e2*weight_smp^0.5 
#w0 <- e0*weight_smp^0.5 
#w1 <- e1*weight_smp^0.5 
#mean(w0^2)
#mean(w2^2)+var(w2)
#mean(w1^2)



transformed_par <- data.frame(e0,weight_smp,framework$smp_data[,framework$smp_domains])
colnames(transformed_par) <- c("e0","weight_smp",framework$smp_domains)
#transformed_par <- data.frame(dep_var,weights_tmp=weight_smp,framework$smp_data)
#transformed_par$ing_lab_pc_v2 <- transformed_par$ing_lab_pc_v2 - indep_smp %*% betas + betas[1]       
random_arg <- NULL 

if (!is.null(framework$smp_subdomains) && !is.null(framework$pop_subdomains)) {
# Do two fold model 
random_arg <- list(as.formula(~1),as.formula(~1))
names(random_arg) <- c(framework$smp_domains,framework$smp_subdomains)
}
else {
random_arg[framework$smp_domains] <- list(as.formula(~1))
}
args <- list(fixed=e0~1,
             data = transformed_par,
             random = random_arg,
             method = framework$nlme_method,weights=~1/weight_smp,
             control = nlme::lmeControl(maxIter = framework$nlme_maxiter,
                                        tolerance = framework$nlme_tolerance,
                                        opt = framework$nlme_opt,
                                        optimMethod = framework$nlme_optimmethod, 
                                        msMaxIter=framework$nlme_msmaxiter,
                                        msTol=framework$nlme_mstol,
                                        returnObject = framework$nlme_returnobject 
             ))
revised_var <- do.call(nlme:::lme,args)

if (!is.null(framework$smp_subdomains) && !is.null(framework$pop_subdomains)) {
  sigmau2est <- as.numeric(VarCorr(revised_var)[2,1])
  dof_adj_h <- (framework$N_subdom_smp-1)/(framework$N_subdom_smp-ncol(indep_smp))
  sigmah2est <- as.numeric(VarCorr(revised_var)[4,1])*dof_adj_h
} 
else {
  sigmau2est <- as.numeric(nlme::VarCorr(revised_var)[1, 1])  
  sigmah2est <- 0 
}

dof_adj_e <- (framework$N_smp-1)/(framework$N_smp-ncol(indep_smp))
sigmae2est <- revised_var$sigma^2*dof_adj_e
dof_adj_u <- (framework$N_dom_smp-1)/(framework$N_dom_smp-ncol(indep_smp))

sigmau2est <- sigmau2est*dof_adj_u  
return(list(sigmae2est=sigmae2est,sigmau2est=sigmau2est))
}



#analytic functions 
# This function calculates the expected value, intially in the case of no transformation 
analytic <- function(transformation,
                     framework,
                     lambda,  
                     shift, 
                     model_par, 
                     gen_model, 
                     fixed, 
                     Ydump
                     ) {

  if(!is.null(framework$aggregate_to_vec)){
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
  
  
  


  
  

# construct vector for variance of random effect, copied from errors_gen line 568
sigma2vu <- vector(length = framework$N_pop)
# variance of random effect for out-of-sample domains
sigma2vu[!framework$obs_dom] <- model_par$sigmau2est
# variance of random effect for in-sample domains
sigma2vu[framework$obs_dom] <- rep(gen_model$sigmav2est,framework$n_pop[framework$dist_obs_dom])
#sigma2vu[framework$obs_dom] <- rep(gen_model$sigmav2est,framework$n_pop[framework$dist_obs_dom])

if (model_par$sigmah2est==0) {
  sigma2eta <-0 
} 
  else {
    sigma2eta <- vector(length=framework$N_pop)
    # variance of subarea random effect for out-of-sample domains    
    sigma2eta[!framework$obs_subdom] <- model_par$sigmah2est 
    sigma2eta[framework$obs_subdom] <- rep(gen_model$sigmai2est,framework$n_pop_subdom[framework$dist_obs_subdom])
  }


var <- sigma2vu+sigma2eta+model_par$sigmae2est



# do mean with function 
indicators <- matrix(ncol=length(framework$indicator_names),nrow=framework$N_pop)
colnames(indicators) <- framework$indicator_names

if ("Mean" %in% framework$indicator_names) {
indicators[,"Mean"] <- expected_untransformed_mean(gen_model$mu,var=var, transformation=transformation,lambda=lambda) 
}
if ("Head_Count" %in% framework$indicator_names) {
indicators[,"Head_Count"] <- expected_head_count(mu=gen_model$mu,var=var, transformation=transformation,lambda=lambda,threshold=framework$threshold,shift=shift)
}
if ("Poverty_Gap" %in% framework$indicator_names) {
if ("Head_Count" %in% framework$indicator_names) {
  Head_Count_temp <- indicators[,"Head_Count"]
  } 
  else {
  Head_Count_temp <- matrix(ncol=1,nrow=framework$N_pop)
  Head_Count_temp <- expected_head_count(mu=gen_model$mu,var=var, transformation=transformation,lambda=lambda,threshold=framework$threshold,shift=shift)
}
  conditional_mean <- conditional_untransformed_mean(mu=gen_model$mu,var=var, transformation=transformation,lambda=lambda,threshold=framework$threshold) 
  indicators[,"Poverty_gap"]<- Head_Count_temp*(1-conditional_mean*Head_Count_temp)/framework$threshold
  }
  

#indicators[,"Median"] <- _percentile(mu=gen_model$mu,var=var, transformation=transformation,lambda=lambda,shift=shift,p=0.5)
#indicators[,"Quantile_10"] <- _percentile(mu=gen_model$mu,var=var, transformation=transformation,lambda=lambda,shift=shift,p=0.1) 
#indicators[,"Quantile_25"] <- _percentile(mu=gen_model$mu,var=var, transformation=transformation,lambda=lambda,shift=shift,p=0.25) 
#indicators[,"Quantile_75"] <- _percentile(mu=gen_model$mu,var=var, transformation=transformation,lambda=lambda,shift=shift,p=0.75) 
#indicators[,"Quantile_90"] <- _percentile(mu=gen_model$mu,var=var, transformation=transformation,lambda=lambda,shift=shift,p=0.9) 


#Output Ydump if selected 
if (!is.null(Ydump)) {
  rand_eff_pop <- rep(model_par$rand_eff, framework$n_pop)
  Ydumpdf <- data.frame(framework$pop_domains_vec,indicators,gen_model$mu,var,sigma2vu,sigma2eta,model_par$sigmae2est,framework$pop_data[,framework$pop_weights],framework$obs_dom,rand_eff_pop,gen_model$mu_fixed)
  colnames(Ydumpdf) <- c("Domain",framework$indicator_names,"Mu","Variance","Var_vu","Var_eta","Var_eps","pop_weight","Observed_dom","Area_reffect","Xbhat")
  write.table(Ydumpdf,file=Ydump,row.names = FALSE,append=FALSE,col.names=T, sep=",")
}






point_estimates <- aggregate_weighted_mean(indicators,by=list("Domain" = pop_domains_vec_tmp),w=pop_weights_vec) 
#point_estimates <- cbind(point_estimates, data.frame(matrix(ncol=9-length(framework$indicator_names),nrow=N_dom_pop_tmp)))
colnames(point_estimates) <- c("Domain",framework$indicator_names)
return(point_estimates)
} # end analytic 

expected_untransformed_mean <- function(mu=mu,var=var,transformation=transformation,lambda=lambda) {
  if (transformation=="no") {
    expected_mean <- mu   
  }
  else if (transformation=="arcsin") {
    term1 <- arcsin_transform_back(mu) # y = sin(x)^2,  
    # dy/dx = 2(sin(x)*cos(x))  
    dy2dx <- -2*sin(mu)^2+2*cos(mu)^2 # product rule for differentiation, f=2sin(x),g=cos(x),f'=2cos(x),g'=-sin(x),dy2dx=f*g'+g*f'=-2sin(x)^2+2*cos(X)^2   
    #dy3dx = -4*sin(mu)*cos(m)-4*cos(x)*sin(x)=-8*sin(x)*cos(x)=-4*dy/dyx
    #dy4/dx <- -4*dy2dx, dy5/dx=-4*dy3dx = 32*sin(x)*cos(x)=16*dydx, dy6/dx <- 16*dy2dx, dy8/dx <- -64*dy2dx
    expected_mean <- term1+0.5*dy2dx*(var)-(1/24)*-4*dy2dx*3*(var^2)+(1/720)*-16*dy2dx*(var^3)+(1/40320)*-64*dy2dx*(var^4)  # 9th order Taylor series
  }
  else if (transformation=="log") {
    expected_mean <- exp(mu+(0.5*var))
  }
  else if (transformation=="log.shift") {
  expected_mean <- exp(mu+(0.5*var))-lambda
  }
  
    return(expected_mean)
}

conditional_untransformed_mean <- function(mu=mu,var=var,transformation=transformation,lambda=lambda,threshold=threshold) {
# first get conditional mean in transformed matric
conditional_transformed_mean <- etruncnorm(a=-Inf,b=threshold-mu,mean=mu,sd=sqrt(var))
conditional_untransformed_mean <- expected_untransformed_mean(mu=conditional_transformed_mean,var=var, transformation=transformation,lambda=lambda) 
return(conditional_untransformed_mean)
}

expected_head_count <- function(mu=mu,threshold=threshold,var=var,transformation=transformation,lambda=lambda,shift=shift) {
  # do poverty by transforming threshold, then applying normal CDF. 
  transformed_threshold <- transformation(y=threshold,transformation=transformation,lambda=lambda,shift=shift)
  expected_head_count <- pnorm(transformed_threshold$y - mu,sd=var^0.5) # formula for head count 
  return(expected_head_count)
}


transformed_percentile <- function(mu=mu,threshold=threshold,var=var,transformation=transformation,lambda=lambda,shift=shift,p=p) {
   mu_percentile <- mu+qnorm(p,0,var)
   transformed_percentile <- back_transformation(y=mu_percentile,transformation=transformation,lambda=lambda,shift=shift) 
   return(transformed_percentile) 
}

# Monte-Carlo approximation ----------------------------------------------------

# The function approximates the expected value (Molina and Rao (2010)
# p.372 (6)). For description of monte-carlo simulation see Molina and
# Rao (2010) p. 373 (13) and p. 374-375
monte_carlo <- function(mixed_model, 
                        transformation,
                        transformation_par, 
                        L,
                        framework,
                        lambda,
                        shift,
                        model_par,
                        gen_model,
                        fixed,
                        Ydump) {

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
  
  for (l in seq_len(L)) {
  
    if (framework$model_parameters=="variable") {
      # variable parameters means they must be drawn every replication
      # so we redraw the parameters  
      # draw error terms 
      # This part is commented out 
      #R <- chol(model_par$varErr) 
      #sigma2<- c(-1,-1)
      #while (sigma2[2]<0 | sigma2[1]<0) {
      #sigma2 <- t(R)  %*% matrix(rnorm(ncol(R)), ncol(R)) + diag(model_par$varErr)
      #model_par$sigmae2est <- sigma2[2]
      #model_par$sigmau2est <- sigma2[1]
      #}
      # draw betas 
      R <- chol(model_par$varFix)
      
      model_par$betas <- t(R)  %*% matrix(rnorm(ncol(R)), ncol(R)) + betas       



      gen_model <- gen_model(
        model_par = model_par,
        fixed = fixed,
        framework = framework,
        dep_var = dep_var 
      )
    } # close condition to redraw parameters 
    
    
    # Errors in generating model: individual error term and random effect
    # See below for function errors_gen.
    errors <- errors_gen(
      framework = framework,
      model_par = model_par,
      gen_model = gen_model
    )

    # Prediction of population vector y
    # See below for function prediction_y.
    population_vector <- prediction_y(
      transformation = transformation,
      lambda = lambda,
      shift = shift,
      errors_gen = errors,
      gen_model = gen_model, 
      framework = framework,
      fixed = fixed
    )

    if(!is.null(framework$pop_weights)){
      pop_weights_vec <- framework$pop_data[[framework$pop_weights]]
    }else{
      pop_weights_vec <- rep(1, nrow(framework$pop_data))
    }
    if (!is.null(Ydump)){
      if (!is.null(framework$smp_subdomains) && !is.null(framework$pop_subdomains)) {
        Ydumpdf <- data.frame(rep(l,nrow(framework$pop_data)), framework$pop_domains_vec,population_vector,gen_model$mu,errors$vu,errors$eta,errors$epsilon)
      }
      else {
        Ydumpdf <- data.frame(rep(l,nrow(framework$pop_data)), framework$pop_domains_vec,population_vector,gen_model$mu,errors$vu,errors$epsilon)  
      }
      
      
      
      
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
  } # End for loop


  # Point estimations of indicators by taking the mean

  point_estimates <- data.frame(
    Domain = unique(pop_domains_vec_tmp),
    apply(ests_mcmc, c(3), rowMeans)
  )
  colnames(point_estimates) <- c("Domain", framework$indicator_names)
  return(point_estimates)
} # End Monte-Carlo


monte_carlo_dt <- function(mixed_model, 
                            transformation,
                            transformation_par, 
                            L,
                            framework,
                            lambda,
                            shift,
                            model_par,
                            gen_model,
                            fixed,
                            Ydump) {

  # Preparing matrices for indicators for the Monte-Carlo simulation
  
  if(!is.null(framework$aggregate_to_vec)){
    N_dom_pop_tmp <- framework$N_dom_pop_agg
    pop_domains_vec_tmp <- framework$aggregate_to_vec
  } else {
    N_dom_pop_tmp <- framework$N_dom_pop
    pop_domains_vec_tmp <- framework$pop_domains_vec
  }
  
  if (framework$model_parameters=="variable") {
    # variable parameters means they must be drawn every replication
    # so we redraw the parameters  
    # draw error terms 
    R <- chol(model_par$varErr) 
    sigma2<- c(-1,-1)
    while (sigma2[2]<0 | sigma2[1]<0) {
      sigma2 <- t(R)  %*% matrix(rnorm(ncol(R)), ncol(R)) + diag(model_par$varErr)
      model_par$sigmae2est <- sigma2[2]
      model_par$sigmau2est <- sigma2[1]
    }
    # draw betas 
    R <- chol(model_par$varFix)
    model_par$betas <- t(R)  %*% matrix(rnorm(ncol(R)), ncol(R)) + model_par$betas      
    
    
    
    gen_model <- gen_model(
      model_par = model_par,
      fixed = fixed,
      framework = framework,
      dep_var = dep_var 
    )
  } # close condition to redraw parameters
  
   errors_mat <- errors_gen_mat(
    framework = framework,
    model_par = model_par, 
    gen_model = gen_model,
    L=L 
  )
  
   population_vector_dt <- prediction_y_dt(
     transformation = transformation,
     lambda = lambda,
     shift = shift,
     errors_gen = errors_mat,
     gen_model = gen_model, 
     framework = framework,
     fixed = fixed,
     L=L 
   )

     
    if (!is.null(Ydump)){
      vu_vec <- matrix(errors_mat$vu,ncol=1)
      epsilon_vec <- matrix(errors_mat$epsilon,ncol=1)
      population_vector_vec <- matrix(population_vector_dt,ncol=1)
      Lindex <-rep(1:L,each=nrow(framework$pop_data))  
      

      
      
      Ydumpdf <- data.frame(Lindex, rep(framework$pop_domains_vec,L),population_vector_vec,rep(gen_model$mu,L),vu_vec,epsilon_vec) 
      colnames(Ydumpdf) <- c("L","Domain","Simulated_Y","XBetahat","eta","epsilon")
      write.table(Ydumpdf,file=Ydump,row.names = FALSE,append=FALSE,col.names=T, sep=",") 
    }
     
   
   
   

       
     y <- cbind(Domain=as.data.table(framework$pop_domain),population_vector_dt) 
     y <- cbind(y,pop_weights=as.data.table(framework$pop_data[[framework$pop_weights]]))
indicators <- function(y,pop_weights,threshold, framework) {
  lapply(framework$indicator_list, function(f) f(y,pop_weights=pop_weights,threshold=framework$threshold))
}

ests_mcmc <- y[,unlist(lapply(.SD, indicators,pop_weights=pop_weights.V1,framework=framework),recursive=FALSE),by=Domain.V1,.SDcols=-ncol(y)]
if (!is.null(framework$indicator_list[["Quantiles"]])) {
  #reshape wide 5 quantiles to columns 
  #add a vector with the quantiles 
    dcast(ests_mcmc,Domain.V1 ~ Quantile_names,value_var = pattern("Quantiles+") )
}



ests_mcmc2 <- data.table::melt(data=ests_mcmc,id.vars = "Domain.V1",measure.vars = patterns(paste0("\\.",framework$indicator_names,"+")))
point_estimates <- data.frame(ests_mcmc2[,lapply(.SD,mean),by=Domain.V1,.SDcols=3:ncol(ests_mcmc2)])
colnames(point_estimates) <- c("Domain", framework$indicator_names)
  return(point_estimates)
} # End Monte-Carlo_vec

errors_gen_mat <- function(framework, model_par, gen_model,L) {
  epsilon <- matrix(rnorm(framework$N_pop*L,0,sqrt(model_par$sigmae2est)),ncol=L)
  
  
  #epsilon <- rnorm(framework$N_pop*L, 0, sqrt(model_par$sigmae2est))
  
  # empty matrix for new random effect in generating model
  vu <- matrix(nrow = framework$N_pop,ncol=L)
  # matrix of random effects  for out-of-sample domains
  vu[!framework$obs_dom,] <- matrix(rep(rnorm(framework$N_dom_unobs*L,0,sqrt(model_par$sigmau2est)),rep(framework$n_pop[!framework$dist_obs_dom],L)),ncol=L)
 
  vu[framework$obs_dom,] <- matrix(rep(rnorm(framework$N_dom_smp*L,0,sqrt(gen_model$sigmav2est)),rep(framework$n_pop[ framework$dist_obs_dom],L)),ncol=L)
  

  #framework$obs_dom_vec <- rep(framework$obs_dom,L)
  #framework$n_pop_vec <- rep(framework$n_pop[!framework$dist_obs_dom],L)
  # new random effect for out-of-sample domains
  #N_dom_unobs_vec <- framework$N_dom_unobs*L
  #vu[!framework$obs_dom_vec] <- rep(
  #  rnorm(
  #    framework$N_dom_unobs*L,
  #    0,
  #    sqrt(model_par$sigmau2est)
  #  ),
  #  framework$n_pop_vec
  #)
  # new random effect for in-sample-domains
  #framework$n_pop_vec <- rep(framework$n_pop[framework$dist_obs_dom],L)
  #vu[framework$obs_dom_vec] <- rep(
  #  rnorm(
  #    framework$N_dom_smp*L,
  #    0,
  #    sqrt(gen_model$sigmav2est)
  #  ),
  #  framework$n_pop_vec
  #)
  # individual error term in generating model epsilon
  
  return(list(epsilon = epsilon, vu = vu))
} # End errors_gen_mat 


prediction_y_dt <- function(transformation,
                         lambda,
                         shift,
                         errors_gen,
                         gen_model, 
                         framework,
                         fixed,
                         L=L) {
  
  
  
  # predicted population income matrix 
  
  mu <- data.table::setDT(rep(list(as.vector(gen_model$mu)),L))
  epsilon <- data.table::as.data.table(errors_gen$epsilon)
  vu <- data.table::as.data.table(errors_gen$vu)
  y_pred <- vu+epsilon+mu 
  

  # back-transformation of predicted population income matrix, by column 
  y_pred [ ,1:ncol(y_pred) := y_pred[,lapply(.SD,back_transformation,transformation = transformation,lambda = lambda,shift = shift,framework = framework,fixed = fixed)]]
  #y_pred <- apply(y_pred,2,back_transformation,transformation = transformation,lambda = lambda,shift = shift,framework = framework,fixed = fixed)
  
  
  #y_pred <- back_transformation(
  #  y = y_pred,
  #  transformation = transformation,
  #  lambda = lambda,
  #  shift = shift,
  #  framework = framework,
  #  fixed = fixed
  #)
  y_pred[mapply(is.infinite, y_pred)] <- NA
  
  #y_pred[!is.finite(y_pred)] <- 0
  return(y_pred)
} # End prediction_y

# The function errors_gen returns error terms of the generating model.
# See Molina and Rao (2010) p. 375 (20)

errors_gen <- function(framework, model_par, gen_model) {

    epsilon <- rnorm(framework$N_pop, 0, sqrt(model_par$sigmae2est))

  # empty vector for new random effect in generating model
  vu <- vector(length = framework$N_pop)
  # new random effect for out-of-sample domains
  vu[!framework$obs_dom] <- rep(
    rnorm(
      framework$N_dom_unobs,
      0,
      sqrt(model_par$sigmau2est)
    ),
    framework$n_pop[!framework$dist_obs_dom]
  )
  # new random effect for in-sample-domains
  vu[framework$obs_dom] <- rep(
    rnorm(
      rep(1, framework$N_dom_smp),
      0,
      sqrt(gen_model$sigmav2est)
    ),
    framework$n_pop[framework$dist_obs_dom]
  )
  # individual error term in generating model epsilon

  return(list(epsilon = epsilon, vu = vu))
} # End errors_gen

# The function prediction_y returns a predicted income vector which can be used
# to calculate indicators. Note that a whole income vector is predicted without
# distinction between in- and out-of-sample domains.
prediction_y <- function(transformation,
                         lambda,
                         shift,
                         errors_gen,
                         gen_model, 
                         framework,
                         fixed) {

  
  
  # predicted population income vector
  y_pred <- gen_model$mu + errors_gen$epsilon + errors_gen$vu

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
} # End prediction_y
