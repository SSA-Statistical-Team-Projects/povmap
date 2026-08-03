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

  call_xgb <- object$out_call

  total_dom <- object$framework$domains_total
  in_dom <- object$framework$domains_in
  oos_dom <- object$framework$domains_out

  dom_info <- data.frame(in_dom, oos_dom, total_dom)
  rownames(dom_info) <- c("")
  colnames(dom_info) <- c("In-sample", "Out-of-sample", "Total")

  smp_size <- object$framework$N_smp
  pop_size <- object$framework$N_pop

  smp_size_dom <- summary(as.data.frame(object$framework$ni_smp)[, "Freq"])
  pop_size_dom <- summary(as.data.frame(object$framework$ni_pop)[, "Freq"])

  sizedom_smp_pop <- rbind(
    Sample_domains = smp_size_dom,
    Population_domains = pop_size_dom
  )





  y <- object$smp_data[,c(object$framework$outcome,object$framework$sub_domains,object$framework$domains,object$framework$smp_weights)]
  #colnames(y)[2] <- "sub_domains"
  y <- na.omit(merge(x=y,y=object$yhat,all.x=T,by=object$framework$sub_domains))
  # Two distinct in-sample fit statistics, reported separately so neither can be
  # mistaken for the other (see print.summary.xgb legend):
  #   squared_correlation = cor(y, yhat)^2, the squared Pearson correlation.
  #     Invariant to a linear rescaling of the prediction; measures co-variation
  #     only. This is the quantity historically labelled "R2".
  #   r2_prop_var = 1 - SSE/SST, the proportion of variance explained. Penalizes
  #     bias and scale compression and can be negative. Numerator and denominator
  #     use matched denominators (both are sums), so it is exactly 1 - SSE/SST.
  squared_correlation <- cor(y[,object$framework$outcome],y$hat)^2
  r2_prop_var <- 1 - sum((y[,object$framework$outcome] - y$hat)^2) /
    sum((y[,object$framework$outcome] - mean(y[,object$framework$outcome]))^2)
  mae <- mean(abs(y[,object$framework$outcome] - y$hat))
  rank_cor <- cor(y[,object$framework$outcome], y$hat, method = "spearman")
  y_yhat <- y[,c(object$framework$outcome,"hat")]
  area_means <- collapse::fmean(x=y_yhat,g=y[,object$framework$domains],w=y[,object$framework$smp_weights])
  area_squared_correlation <- cor(area_means[,object$framework$outcome],area_means$hat)^2
  area_r2_prop_var <- 1 - sum((area_means[,object$framework$outcome] - area_means$hat)^2) /
    sum((area_means[,object$framework$outcome] - mean(area_means[,object$framework$outcome]))^2)
  area_mae <- mean(abs(area_means[,object$framework$outcome] - area_means$hat))
  area_rank_cor <- cor(area_means[,object$framework$outcome], area_means$hat, method = "spearman")

  # Variance decomposition: area effect (sigma2_v) and idiosyncratic (sigma2_e)
  y$resid <- y[, object$framework$outcome] - y$hat
  domain_mean_resid <- collapse::fmean(
    x = y$resid,
    g = y[, object$framework$domains],
    w = y[, object$framework$smp_weights]
  )
  sigma2_v <- var(domain_mean_resid)
  y$domain_mean_resid <- domain_mean_resid[match(y[, object$framework$domains], names(domain_mean_resid))]
  sigma2_e <- var(y$resid - y$domain_mean_resid)

  var_decomp <- data.frame(
    Sigma2_v = sigma2_v,
    Sigma2_e = sigma2_e,
    row.names = ""
  )

  coeff_det <- data.frame(
    Squared_correlation      = squared_correlation,
    R2_prop_var              = r2_prop_var,
    MAE                      = mae,
    Rank_cor                 = rank_cor,
    Area_squared_correlation = area_squared_correlation,
    Area_R2_prop_var         = area_r2_prop_var,
    Area_MAE                 = area_mae,
    Area_Rank_cor            = area_rank_cor,
    row.names                = ""
  )

  feature_names <- xgboost:::xgb.feature_names(object$model)
  n_features <- length(feature_names)

  # information on xgb:
  xgb_info <- data.frame(c(
    object$transformation,
    xgboost::xgb.get.num.boosted.rounds(object$model),
    attributes(object$model)$params$max_depth,
    n_features)
  )

  colnames(xgb_info) <- NULL
  rownames(xgb_info) <- c(
    "Transformation","Number of boosting interations:", "Maximum depth of a tree:",
    "Number of independent variables:"
  )




  sum_xgb <- list(
    call_xgb = call_xgb,
    dom_info = dom_info,
    smp_size = smp_size,
    pop_size = pop_size,
    sizedom_smp_pop = sizedom_smp_pop,
    coeff_determ = coeff_det,
    var_decomp = var_decomp,
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
  cat("Fit statistics (in-sample):\n")
  cat("________________________________________________________________\n")
  cat("Squared_correlation: squared Pearson correlation, cor(y, yhat)^2\n")
  cat("R2_prop_var: proportion of variance explained, 1 - SSE/SST (can be negative)\n")
  cat("Prefix Area_: computed on the domain (area) means rather than the units\n")
  print(x$coeff_determ)
  cat("\n")
  cat("Variance components:\n")
  cat("________________________________________________________________\n")
  cat("Sigma2_v: area effect variance\n")
  cat("Sigma2_e: idiosyncratic error variance\n")
  print(x$var_decomp)
  cat("\n")
  cat("Boosting component: \n")
  cat("________________________________________________________________\n")
  print(x$xgb_info)
  cat("\n")
}

