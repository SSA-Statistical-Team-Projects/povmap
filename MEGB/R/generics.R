
#' Coef Method for megb Class
#' @param megb_object object of class MEGB containing the boosting model
#' @param ... additional arguments
#' @return List containing random effects coefficients and feature importance
#' @export
coef.megb <- function(megb_object, ...) {
  UseMethod("coef", megb_object)
}


#' Default Coef Method for MEGB Class
#' @param megb_object object of class MEGB containing the boosting model
#' @param ... additional arguments
#' @return List of random effects coefficients and feature importance
#' @export
coef.MEGB <- function(megb_object, ...) {
  list(
    feature_importance = megb_object$megb_model$importance_matrix,
    random_effects = unclass(coef(megb_object$megb_model$effect_model))[[1]]
  )
}

#' Save objects of class MEGB
#' @param megb_object object of class MEGB containing the boosting model
#' @param file file to save object to
#' @param ... additional arguments
#' @return NULL invisibly
#' @export
save_megb <- function(megb_object, file, ...) {
  saveRDS(megb_object, file = file)
}

#' Load objects of class MEGB
#' @param file file to load object from
#' @param ... additional arguments
#' @return an R object
#' @export
load_megb <- function(file, ...) {
readRDS(file)
}

#' S3 print-method for MEGB-objects
#' @param megb_object object of class MEGB containing the boosting model
#' @param ... additional arguments
#' @return summary of model information
#' @export

print.MEGB <- function(megb_object, ...) {
  cat("________________________________________________________________
")
  cat("Mixed Effects Gradient Boosting (MEGB) Model
")
  cat("________________________________________________________________
")
  cat("
")

  # Model Call
  cat("Model Call:
")
  print(megb_object$call)
  cat("
")

  # Summary of most important model information
  cat("Model Summary:
")
  cat(sprintf(
    "Number of Covariates (Fixed Effects): %d |\n Number of Random Effects (Domains): %d |\n Number of EM Iterations Used: %d
",
    length(megb_object$megb_model$boosting$feature_names),
    nrow(megb_object$Indicators),
    megb_object$megb_model$iterations_used
  ))
  cat("
")

  # Unit-level Diagnosis
  cat("Unit-Level Residual Diagnostics:
")
  cat("Residual Standard Deviation: ", round(sd(megb_object$megb_model$residuals), 3), "
")
  cat("Residual Summary:
")
  print(summary(megb_object$megb_model$residuals))
  cat("
")

  # RMSE
  cat("Estimated Root Mean Squared Error (RMSE) Summary:
")
  if(is.null(megb_object$MSE_Estimates)){
    cat("MSE was not calculated.
")
  } else {
    print(summary(sqrt(megb_object$MSE_Estimates$MSE_estimates[, 2])))
  }
  cat("
")

  # Boosting Information
  cat("Boosting Information:
")
  cat(sprintf("Number of Boosting Iterations: %d
", megb_object$megb_model$boosting$niter))

}


