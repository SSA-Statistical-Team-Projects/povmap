# Summarizes an povmap ell Object

#' @export
#' @importFrom moments skewness kurtosis
#' @importFrom MuMIn r.squaredGLMM
#' @rdname povmap_summaries

# importFrom nlme predict.lme -- old code 

summary.ell <- function(object, ...) {
  throw_class_error(object, "ell")

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
  residuals <- object$model_par$residuals-rep(object$model_par$mean_residuals,object$framework$n_smp)
  
  
  
  skewness_res <- skewness(residuals)
  kurtosis_res <- kurtosis(residuals)
  variance_res <- object$model_par$sigmae2est 
  variance_ran <- object$model_par$sigmau2est 
  
  
  
    #ranef <- ranef(object$model)$"(Intercept) # traditional defition, does not account for weights 
    dist_obs_dom <- unique(object$framework$pop_domains_vec) %in% unique(object$framework$smp_domains_vec)
    ranef <- object$model_par$mean_residuals 
    
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
 

  r_squared <- summary(object$model)$r.squared

  groups=object$framework$smp_data[,object$framework$smp_domains]
  if (is.null(object$framework$weights)) {
    temp_weights <- rep(1,object$framework$N_smp)
  }
  else {
    temp_weights <- object$framework$smp_data[,object$framework$weights]
  }
  y_yhat <- data.frame("Y" =  object$model$model[,1],"Yhat" = object$model$model[,1] - object$model_par$residuals, 
                         weights=temp_weights)
  y_yhat_bar <- aggregate_weighted_mean(y_yhat,by=list(groups),w=temp_weights)
  
  r2_area <- cor(y_yhat_bar[,2],y_yhat_bar[,3])^2
  

  
  coeff_det <- data.frame(
    R2    = r_squared[1],
    Area_R2 = r2_area,
    row.names = ""
  )
  

    

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
    coeff_determ = coeff_det,
    call = call_povmap
  )

  class(sum_povmap) <- c("summary.ell", "povmap")
  sum_povmap
}


#' @export
print.summary.ell <- function(x, ...) {
  throw_class_error(x, "ell")
  cat("ELL estimation\n")
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
    cat("Explanatory measures for the random effects model:\n")
  }
  print(x$coeff_determ)
  cat("\n")
  if (is.null(x$call$weights)) {
    cat("Residual diagnostics:\n")
  } else {
    cat("Residual diagnostics for the random effects model:\n")
  }
  print(x$normality)
  cat("\n")
  cat("Estimated variance of random effects:\n")
  print(x$variance)
  cat("\n")
  if (is.null(x$transform)) {
    cat("Transformation: No transformation \n")
  } else {
    cat("Transformation:\n")
    print(x$transform)
  }
}


