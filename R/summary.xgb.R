#' Summarize an extreme gradient boosting model for domain-level averages
#'
#' Additional information about the data, model and components of an \code{xgb} object
#' are extracted. The returned object is suitable for printing with \code{print}.
#'
#'
#' @param object an object of class \code{xgb}, \code{emdi}, containing point
#' estimates, MSE, and confidence interval estimates.
#' @param ... optional additional inputs that are ignored for this method.
#'
#' @return An object of class \code{summary.xgb} including information about the sample
#' and population data and extreme gradient boosting specific metrics.
#'
#' @export
#'
#' @examples
#' \donttest{
#' # Loading data - population and sample data
#' data("eusilcA_pop")
#' data("eusilcA_smp")
#'
#' # Create subdomains; equal to individuals in each area
#' eusilcA_smp$subDomain <- ave(eusilcA_smp$district,
#'                              eusilcA_smp$district,
#'                              FUN = seq_along)
#' eusilcA_pop$subDomain <- ave(eusilcA_pop$district,
#'                              eusilcA_pop$district,
#'                              FUN = seq_along)
#'
#' xgb_model <- xgb(fixed = eqIncome ~ eqsize + cash + self_empl +
#'                  unempl_ben + age_ben + surv_ben + sick_ben + dis_ben +
#'                  rent + fam_allow + house_allow + cap_inv + tax_adj +
#'                  district + subDomain,
#'                  smp_data = eusilcA_smp,
#'                  pop_data = eusilcA_pop,
#'                  domains = "district",
#'                  sub_domains = "subDomain")
#'
#' # Receive first overview
#' summary(xgb_model)
#'}

summary.xgb <- function(object, ...) {

  call_xgb <- object$xgbModel$call

  total_dom <- object$xgbModel$domains_total
  in_dom <- object$xgbModel$domains_in
  oos_dom <- object$xgbModel$domains_out

  dom_info <- data.frame(in_dom, oos_dom, total_dom)
  rownames(dom_info) <- c("")
  colnames(dom_info) <- c("In-sample", "Out-of-sample", "Total")

  smp_size <- object$xgbModel$N_smp
  pop_size <- object$xgbModel$N_pop

  smp_size_dom <- summary(as.data.frame(object$xgbModel$ni_smp)[, "Freq"])
  pop_size_dom <- summary(as.data.frame(object$xgbModel$ni_pop)[, "Freq"])

  sizedom_smp_pop <- rbind(
    Sample_domains = smp_size_dom,
    Population_domains = pop_size_dom
  )


  # information on xgb:
  xgb_info <- data.frame(c(
    object$xgbModel$transformation,
    object$xgbModel$niter,
    object$xgbModel$params$max_depth,
    object$xgbModel$nfeatures)
  )

  colnames(xgb_info) <- NULL
  rownames(xgb_info) <- c(
    "Transformation","Number of booting interations:", "Maximum depth of a tree:",
    "Number of independent variables:"
  )


  sum_xgb <- list(
    call_xgb = call_xgb,
    dom_info = dom_info,
    smp_size = smp_size,
    pop_size = pop_size,
    sizedom_smp_pop = sizedom_smp_pop,
    xgb_info = xgb_info
  )

  class(sum_xgb) <- c("summary.xgb", "emdi")
  sum_xgb
}

# Generic print function for summary.xgb --------------------------------------------
#' @export
print.summary.xgb <- function(x, ...) {
  #class_error(object = x)
  cat("________________________________________________________________\n")
  cat("Extreme Gradient Boosting for Small Area Estimation\n")
  cat("________________________________________________________________\n")
  cat("Call:\n")
  print(x$call_xgb)
  cat("\n")
  cat("Domains\n")
  cat("________________________________________________________________")
  cat("\n")
  print(x$dom_info)
  cat("\n")
  cat("Totals:\n")
  cat("Units in sample:", x$smp_size, "\n")
  if (!is.null(x$pop_size)) {
    cat("Units in population:", x$pop_size, "\n")
  }
  cat("\n")
  print(x$sizedom_smp_pop)
  cat("\n")
  cat("Boosting component: \n")
  cat("________________________________________________________________\n")
  print(x$xgb_info)
  cat("\n")
}

