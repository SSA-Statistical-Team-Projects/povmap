# Summarizes an povmap ebp Object

#' @export
#' @importFrom moments skewness kurtosis
#' @importFrom MuMIn r.squaredGLMM
#' @rdname povmap_summaries

# importFrom nlme predict.lme -- old code 

summary.ebp <- function(object, ...) {
  throw_class_error(object, "ebp")

  call_povmap <- object$call

  N_dom_unobs <- object$framework$N_dom_unobs
  N_dom_smp <- object$framework$N_dom_smp

  N_subdom_unobs <-   object$framework$N_subdom_unobs
  N_subdom_smp <- object$framework$N_subdom_smp
  
  smp_size <- object$framework$N_smp
  pop_size <- object$framework$N_pop

  smp_size_dom <- summary(as.data.frame(
    table(object$framework$smp_domains_vec)
  )[, "Freq"])
  pop_size_dom <- summary(as.data.frame(
    table(object$framework$pop_domains_vec)
  )[, "Freq"])
  sizedom_smp_pop <- rbind(
    Sample_domains = smp_size_dom,
    Population_domains = pop_size_dom
  )

  if (object$transformation == "box.cox" || object$transformation == "dual") {
    transform_method <- data.frame(
      Transformation = object$transformation,
      Method = object$method,
      Optimal_lambda =
        object$transform_param$optimal_lambda,
      Shift_parameter =
        round(object$transform_param$shift_par, 3),
      row.names = ""
    )
  } else if (object$transformation == "log.shift") {
    transform_method <- data.frame(
      Transformation = object$transformation,
      Method = object$method,
      Optimal_lambda =
        object$transform_param$optimal_lambda,
      row.names = ""
    )
  } else if (object$transformation == "log") {
    transform_method <- data.frame(
      Transformation = object$transformation,
      Shift_parameter =
        round(object$transform_param$shift_par, 3),
      row.names = ""
    )
  } else if (object$transformation == "ordernorm") {
    transform_method <- data.frame(Transformation  = object$transformation,
                                   Shift_parameter = 0,
                                   row.names       = ""
    )
  } else if (object$transformation == "arcsin") {
    transform_method <- data.frame(Transformation  = object$transformation,
                                   Shift_parameter = 0,
                                   row.names       = ""
    )
  } else if (object$transformation == "no") {
    transform_method <- NULL
  }
  # traditionally povmap uses this definition of residuals 
  # but it doesn't account for weights properly when using Guadarrama weights or hybrid weights 
  #residuals <- residuals(object$model level = 0,  type = "pearson")
  #raneff <- 
  residuals <- object$model_par$e1
  
  
  
  skewness_res <- skewness(residuals)
  kurtosis_res <- kurtosis(residuals)
  variance_res <- object$model_par$sigmae2est 
  variance_ran <- object$model_par$sigmau2est 
  
  if (!is.null(object$framework$smp_subdomains) && !is.null(object$framework$pop_subdomains))
  { # two fold model
    residuals <- object$model_par$e2
    dist_obs_dom <- unique(object$framework$pop_domains_vec) %in% unique(object$framework$smp_domains_vec)
    dist_obs_subdom <- unique(object$framework$pop_subdomains_vec) %in% unique(object$framework$smp_subdomains_vec)
  ranef <- object$model_par$rand_eff[dist_obs_dom] 
  ranef_sub <- object$model_par$rand_eff_h[dist_obs_subdom] 
   skewness_ran <- skewness(ranef) 
   kurtosis_ran <- kurtosis(ranef)
   skewness_ran_sub <- skewness(ranef_sub)
   kurtosis_ran_sub <- kurtosis(ranef_sub)
   variance_ran_sub <- object$model_par$sigmah2est 
   
   if (length(residuals) > 3 &&
       length(residuals) < 5000) {
     shapiro_p_res <-
       shapiro.test(residuals)[[2]] 
     shapiro_W_res <-
       shapiro.test(residuals)[[1]]
     
   }
   else {
     warning(strwrap(prefix = " ", initial = "",
                     "Number of observations exceeds 5000 or is lower then 3 and
                    thus the Shapiro-Wilk test is not applicable for the
                    residuals."))
     shapiro_W_res <- NA
     shapiro_p_res <- NA
   }
   if (length(ranef) > 3 &&
       length(ranef) < 5000) {
     shapiro_p_ran <- shapiro.test(ranef)[[2]]
     shapiro_W_ran <- shapiro.test(ranef)[[1]]
   } else {
     warning(strwrap(prefix = " ", initial = "",
                     "Number of domains exceeds 5000 or is lower then 3 and thus
                    the Shapiro-Wilk test is not applicable for the area random
                    effects."))
     shapiro_p_ran <- NA
     shapiro_W_ran <- NA
   }
   if (length(ranef_sub) > 3 &&
       length(ranef_sub) < 5000) {
     shapiro_p_ran_sub <- shapiro.test(ranef_sub)[[2]]
     shapiro_W_ran_sub <- shapiro.test(ranef_sub)[[1]]
   } else {
     warning(strwrap(prefix = " ", initial = "",
                     "Number of subdomains exceeds 5000 or is lower then 3 and thus
                    the Shapiro-Wilk test is not applicable for the subarea random
                    effects."))
     shapiro_p_ran_sub <- NA
     shapiro_W_ran_sub <- NA
   }
   
   norm <- data.frame(
     Skewness = c(skewness_res, skewness_ran,skewness_ran_sub), 
     Kurtosis = c(kurtosis_res, kurtosis_ran,kurtosis_ran_sub),  
     Shapiro_W = c(shapiro_W_res, shapiro_W_ran, shapiro_W_ran_sub),
     Shapiro_p = c(shapiro_p_res, shapiro_p_ran, shapiro_p_ran_sub),
     row.names = c("Error", paste0(object$framework$smp_domains," random effect"),paste0(object$framework$smp_subdomains," random effect"))
   )
   var <- data.frame(Variance = c(variance_res,variance_ran,variance_ran_sub),
    row.names = c("Error", paste0(object$framework$smp_domains," random effect"),paste0(object$framework$smp_subdomains," random effect"))                         
   )
  } # close two fold model 
  else {
    #one fold model
    #ranef <- ranef(object$model)$"(Intercept) # traditional defition, does not account for weights 
    dist_obs_dom <- unique(object$framework$pop_domains_vec) %in% unique(object$framework$smp_domains_vec)
    ranef <- object$model_par$rand_eff[dist_obs_dom] 
    
  skewness_ran <- skewness(ranef)
  kurtosis_ran <- kurtosis(ranef)
  if (length(residuals) > 3 &&
    length(residuals) < 5000) {
    shapiro_p_res <-
      shapiro.test(residuals)[[2]]
    shapiro_W_res <-
      shapiro.test(residuals)[[1]]
  } else {
    warning(strwrap(prefix = " ", initial = "",
                    "Number of observations exceeds 5000 or is lower then 3 and
                    thus the Shapiro-Wilk test is not applicable for the
                    residuals."))
    shapiro_p_res <- NA
    shapiro_W_res <- NA
  }

  if (length(ranef) > 3 &&
    length(ranef) < 5000) {
    shapiro_p_ran <- shapiro.test(ranef)[[2]]
    shapiro_W_ran <- shapiro.test(ranef)[[1]]
  } else {
    warning(strwrap(prefix = " ", initial = "",
                    "Number of domains exceeds 5000 or is lower then 3 and thus
                    the Shapiro-Wilk test is not applicable for the random
                    effects."))
    shapiro_p_ran <- NA
    shapiro_W_ran <- NA
  }

  norm <- data.frame(
    Skewness = c(skewness_res, skewness_ran),
    Kurtosis = c(kurtosis_res, kurtosis_ran),
    Shapiro_W = c(shapiro_W_res, shapiro_W_ran),
    Shapiro_p = c(shapiro_p_res, shapiro_p_ran),
    row.names = c("Error", "Random_effect")
  )
  var <- data.frame(Variance = c(variance_res,variance_ran),
                    row.names = c("Error", "Random_effect")
                    )
  } # close one fold model  
  tempMod <- object$model
  tempMod$call$fixed <- object$fixed
  r_squared <- suppressWarnings(r.squaredGLMM(tempMod))
  if (is.matrix(r_squared)) {
    r_marginal <- r_squared[1, 1]
    r_conditional <- r_squared[1, 2]
  } else {
    r_marginal <- r_squared[1]
    r_conditional <- r_squared[2]
  }
  icc_mixed <- icc(object)

  groups=tempMod$data[,object$framework$smp_domains]
  if (is.null(object$framework$weights)) {
    object$framework$weights <- "temp_weights"
    object$model$data$temp_weights <- 1
  }
  y_yhat <- data.frame("Y" = tempMod$data[,1],"marginal" = tempMod$fitted[,1], "conditional"=tempMod$fitted[,2],weights=object$model$data[,object$framework$weights])
  y_yhat_bar <- aggregate_weighted_mean(y_yhat,by=list(groups),w=object$model$data[,object$framework$weights])[,-c(1,5)]
  
  r2_area <- cor(y_yhat_bar[])[2:3,1]^2
  r_marginal_area <- r2_area[1]
  r_conditional_area <- r2_area[2]
  
  
  
  
  coeff_det <- data.frame(
    Marginal_R2    = r_marginal,
    Conditional_R2 = r_conditional,
    Marginal_Area_R2 = r_marginal_area,
    Conditional_Area_R2 = r_conditional_area, 
    row.names      = ""
  )
  
  if (is.null(object$model_par$gamma_sub)) {
    shrinkage <- data.frame(
      gamma = t(as.matrix(summary(object$model_par$gamma))),
      row.names = paste0(object$framework$smp_domains,"")
    )
  }  
  else {
    shrinkage <- data.frame(
      gamma = t(as.matrix(cbind(summary(object$model_par$gamma),summary(object$model_par$gamma_sub)))),
      row.names = c(paste0(object$framework$smp_domains,""),paste0(object$framework$smp_subdomains,""))
      )
  }

    
    
  
  
  

  sum_povmap <- list(
    out_of_smp = N_dom_unobs,
    in_smp = N_dom_smp,
    out_of_smp_sub = N_subdom_unobs, 
    in_smp_sub = N_subdom_smp, 
    size_smp = smp_size,
    size_pop = pop_size,
    size_dom = sizedom_smp_pop,
    smp_size_tab = NULL,
    transform = transform_method,
    normality = norm,
    variance = var, 
    shrinkage = shrinkage,  
    icc = icc_mixed,
    coeff_determ = coeff_det,
    call = call_povmap
  )

  class(sum_povmap) <- c("summary.ebp", "povmap")
  sum_povmap
}


#' @export
print.summary.ebp <- function(x, ...) {
  throw_class_error(x, "ebp")
  cat("Empirical Best Prediction\n")
  cat("\n")
  cat("Call:\n ")
  print(x$call)
  cat("\n")
  cat("Out-of-sample domains: ", x$out_of_smp, "\n")
  cat("In-sample domains: ", x$in_smp, "\n")
  if (!is.null(x$in_smp_sub)) {
  cat("Out-of-sample subdomains: ", x$out_of_smp_sub, "\n")
  cat("In-sample subdomains: ", x$in_smp_sub, "\n")  
  }
  cat("\n")
  cat("Sample sizes:\n")
  cat("Units in sample: ", x$size_smp, "\n")
  cat("Units in population: ", x$size_pop, "\n")
  print(x$size_dom)
  cat("\n")
  if (is.null(x$call$weights)) {
    cat("Explanatory measures:\n")
  } else {
    cat("Explanatory measures for the mixed model:\n")
  }
  print(x$coeff_determ)
  cat("\n")
  if (is.null(x$call$weights)) {
    cat("Residual diagnostics:\n")
  } else {
    cat("Residual diagnostics for the mixed model:\n")
  }
  print(x$normality)
  cat("\n")
  cat("Estimated variance of random effects:\n")
  print(x$variance)
  cat("\n")
  cat("ICC: ", x$icc, "\n")
  cat("\n")
  cat("Shrinkage factors","\n")
  print(x$shrinkage)
  cat("\n")
  if (is.null(x$transform)) {
    cat("Transformation: No transformation \n")
  } else {
    cat("Transformation:\n")
    print(x$transform)
  }
}


#  ICC

icc <- function(object) {
  u <- object$model_par$sigmau2est 
  e <- (object$model_par$sigmae2est+object$model_par$sigmah2est)
  u / (u + e)
}
