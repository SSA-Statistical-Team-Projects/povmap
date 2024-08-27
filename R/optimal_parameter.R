optimal_parameter <- function(generic_opt,
                              fixed,
                              smp_data,
                              smp_domains,
                              transformation,
                              interval,
                              framework) {

if (transformation != "no" && transformation != "log" &&
      transformation != "ordernorm" && transformation != "arcsin"
    && transformation != "logit") {
    # no lambda -> no estimation -> no optmimization

    if (transformation == "box.cox" && any(interval == "default")) {
      interval <- c(-1, 2)
    } else if (transformation == "dual" && any(interval == "default")) {
      interval <- c(0, 2)
    } else if (transformation == "log.shift" && any(interval == "default")) {
      # interval = c(min(smp_data[paste(fixed[2])]),
      # max(smp_data[paste(fixed[2])]))
      span <- range(smp_data[paste(fixed[2])])
      if ((span[1] + 1) <= 1) {
        lower <- abs(span[1]) + 1
      } else {
        lower <- 0
      }

      upper <- diff(span) / 2

      interval <- c(lower, upper)
    }

    # Estimation of optimal lambda parameters

  if (!is.null(framework$nlme_method)) {
    #EBP 
      optimal_parameter <- optimize(generic_opt,
      fixed          = fixed,
      smp_data       = smp_data,
      smp_domains    = smp_domains,
      transformation = transformation,
      interval       = interval,
      framework      = framework,
      maximum        = FALSE
    )$minimum
  } else {
    #ELL 
  optimal_parameter <- optimize(generic_opt_ell,
      fixed          = fixed,
      smp_data       = smp_data,
      smp_domains    = smp_domains,
      transformation = transformation,
      interval       = interval,
      framework      = framework,
      maximum        = FALSE
  )$minimum
        } 
  }
  else {
    optimal_parameter <- NULL
  }

  return(optimal_parameter)
} # End optimal parameter


# Internal documentation -------------------------------------------------------

# Function generic_opt provides estimation method reml to specifiy
# the optimal parameter lambda. Here its important that lambda is the
# first argument because generic_opt is given to optimize. Otherwise,
# lambda is missing without default.

generic_opt <- function(lambda,
                        fixed,
                        smp_data,
                        smp_domains,
                        transformation,
                        framework) {

  # Definition of optimization function for finding the optimal lambda
  # Preperation to easily implement further methods here
  optimization <- if (TRUE) {
    reml(
      fixed = fixed,
      smp_data = smp_data,
      smp_domains = smp_domains,
      transformation = transformation,
      lambda = lambda,
      framework = framework
    )
  }
  return(optimization)
}

generic_opt_ell <- function(lambda,
                        fixed,
                        smp_data,
                        smp_domains,
                        transformation,
                        framework) {
  
  # Definition of optimization function for finding the optimal lambda
  # Preperation to easily implement further methods here
  optimization <- if (TRUE) {
    ml_plm(
      fixed = fixed,
      smp_data = smp_data,
      smp_domains = smp_domains,
      transformation = transformation,
      lambda = lambda,
      framework = framework
    )
  }
  return(optimization)
}



# REML method ------------------------------------------------------------------

reml <- function(fixed = fixed,
                 smp_data = smp_data,
                 smp_domains = smp_domains,
                 transformation = transformation,
                 lambda = lambda,
                 framework = framework) {

  sd_transformed_data <- std_data_transformation(
    fixed = fixed,
    smp_data = smp_data,
    transformation =
      transformation,
    lambda = lambda
  )

  model_REML <- NULL
  try(
    if(!is.null(framework$weights) && framework$weights_type == "nlme_lambda") {

      sd_transformed_data$weights_scaled <- smp_data[,framework$weights] /
        mean(smp_data[,framework$weights], na.rm = TRUE)

      model_REML <- lme(
        fixed = fixed,
        data = sd_transformed_data,
        random = as.formula(paste0("~ 1 | as.factor(", smp_domains, ")")),
        method = framework$nlme_method,
        control = nlme::lmeControl(maxIter = framework$nlme_maxiter,
                                   tolerance = framework$nlme_tolerance,
                                   opt = framework$nlme_opt,
                                   optimMethod = framework$nlme_optimmethod, 
                                   msMaxIter=framework$nlme_msmaxiter,
                                   msTol=framework$nlme_mstol,
                                   returnObject = framework$nlme_returnobject 
        ),
        keep.data = FALSE,
        weights =
          varComb(varIdent(as.formula(
                      paste0("~ 1 | as.factor(", framework$smp_domains, ")"),)),
                  varFixed(as.formula(paste0("~1/", "weights_scaled")))
        )
      )
    } else {
      model_REML <- lme(
        fixed = fixed,
        data = sd_transformed_data,
        random = as.formula(paste0("~ 1 | as.factor(", smp_domains, ")")),
        method = framework$nlme_method,
        control = nlme::lmeControl(maxIter = framework$nlme_maxiter,
                                   tolerance = framework$nlme_tolerance,
                                   opt = framework$nlme_opt,
                                   optimMethod = framework$nlme_optimmethod, 
                                   msMaxIter=framework$nlme_msmaxiter,
                                   msTol=framework$nlme_mstol,
                                   returnObject = framework$nlme_returnobject 
        ),
        keep.data = FALSE
      )
    }, silent = TRUE)

  if (is.null(model_REML)) {
    stop(strwrap(prefix = " ", initial = "",
                 "The likelihood did not converge when estimating an nlme model to select the optimal
                 transformation parameter. Try adjusting nlme_opt, nlme_maxiter, nlme_tolerance, 
                 nlme_optimmethod, or nlme_returnobject when calling ebp, or using a non-adapative transformation. See also
                 help(ebp)."))
  } else {
    model_REML <- model_REML
  }


  log_likelihood <- -logLik(model_REML)

  return(log_likelihood)
}

ml_plm <- function(fixed = fixed,
                 smp_data = smp_data,
                 smp_domains = smp_domains,
                 transformation = transformation,
                 lambda = lambda,
                 framework = framework) {
  
  sd_transformed_data <- std_data_transformation(
    fixed = fixed,
    smp_data = smp_data,
    transformation =
      transformation,
    lambda = lambda
  )
  
  model_PLM <- NULL
  weights_arg <- NULL
  if (!is.null(framework$weight)) {
    weights_arg <- framework$smp_data[,framework$weights]
}

  args <- list(formula=fixed, 
               data = sd_transformed_data, 
               weights = weights_arg ,
               model="random",
               index = framework$smp_domains,
               random.method=framework$random_method)
  
  try(
    model_PLM <- do.call(plm:::plm, args)
, silent = TRUE)

  if (is.null(model_PLM)) {
    stop(strwrap(prefix = " ", initial = "",
                 "The likelihood did not converge when estimating a random effects model to select the optimal
                 transformation parameter. Try using a non-adapative transformation."))
  } else {
    model_PLM <- model_PLM
  }
  log_likelihood <- logLik.plm(model_PLM)[1]
  return(log_likelihood)
}

  # taken from Stata XT manual help for xtreg, re
  logLik.plm <- function(object){
    e=object$residuals
    e2=object$residuals^2
    tosum <- cbind(e2,rep(1,nobs(object)),e)
    group <- plm:::index.panelmodel(object)
    sum <- aggregate(tosum,by=list(group[,1]),FUN=sum) 
    s2e <- object$ercomp$sigma2[1]
    s2u <- object$ercomp$sigma2[2]
    Ti <- sum[,3]
    term1 <- (1/s2e)*(sum[,2]-s2u/(Ti*s2u+s2e)*sum[,4]^2)
    term2 <- log(Ti*s2u/s2e+1)+Ti*log(2*pi*s2e)
    ll <- sum(-0.5*(term1+term2))
    attr(ll,"df") <- nobs(object) - object$df.residual
    attr(ll,"nobs") <- plm::nobs(object)
    return(ll)
  }
  