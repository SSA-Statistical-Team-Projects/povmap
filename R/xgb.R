#' Extreme gradient boosting for domain-level averages
#'
#' The function \code{xgb} employs the extreme gradient boosting methodology introduced
#' by \cite{Merfeld and Newhouse (2023)} to estimate domain-level averages, particularly
#' for small area estimation (SAE) applications. Moreover, to estimate the mean squared
#' error (MSE), a nonparametric residual bootstrap approach is utilized, as described
#' in \cite{Krennmair and Schmid (2022)} and \cite{Merfeld and Newhouse (2023)}.
#'
#' @param fixed a two-sided linear formula object describing the
#' fixed-effects part of the model with the dependent variable on the left
#' of a ~ operator and the explanatory variables on the right, separated
#' by + operators. All variables (except for \code{domains} and \code{subdomains})
#' must be numeric.
#' @param smp_data a data frame that needs to comprise all variables including
#' \code{domains} and \code{sub_domains}.
#' @param smp_weights a character string containing the name of the variable that
#' indicates weights in the \code{smp_data}. The variable has to be numeric.
#' Defaults to \code{NULL}.
#' @param pop_data a data frame that needs to comprise all variables including
#' \code{domains} and \code{sub_domains}.
#' @param pop_weights a character string containing the name of the variable that
#' indicates population weights in \code{pop_data}. The variable has to be
#'  numeric. Defaults to \code{NULL}.
#' @param domains a character string containing the name of a variable
#' that indicates domains in \code{smp_data} and \code{pop_data}. The variable can be
#' numeric or a factor.
#' @param sub_domains character string specifying the variable name that denotes
#' sub-domains within the dataset. This variable must have unique values across
#' observations.
#' @param transformation a character string. Two different transformation
#' types for the dependent variable can be chosen (i) no transformation ("no");
#' (ii) log transformation ("log"); (iii) Arcsin transformation ("arcsin").
#' Defaults to \code{"no"}.
#' @param B a number determining the number of bootstrap populations in the
#' nonparametric residual bootstrap approach used in the MSE estimation. The
#' number must be greater than 1. Defaults to 100. For practical applications,
#' values larger than 200 are recommended.
#' @param conf_level confidence level for the confidence interval. Defaults to 0.95.
#' @param nround maximum number of boosting iterations. Defaults to 100.
#' @param max_depth maximum depth of a tree. Increasing this value will result
#' in a more complex model, increasing the likelihood of overfitting. A value of
#'  0 indicates no limit on the depth. Defaults to 4.
#' @param colsample_bytree subsample ratio of columns when constructing
#' each tree. Subsampling occurs once for every tree constructed. Range of (0, 1].
#' Defaults to 0.6.
#' @param colsample_bylevel subsample ratio of columns for each level.
#' Subsampling occurs once for every new depth level reached in a tree. Columns
#' are subsampled from the set of columns chosen for the current tree. Range of (0, 1].
#' Defaults to 0.6.
#' @param colsample_bynode subsample ratio of columns for each node (split).
#' Subsampling occurs once every time a new split is evaluated. Columns are
#' subsampled from the set of columns chosen for the current level. Range of (0, 1].
#' Defaults to 0.6.
#' @param subsample subsample ratio of the training instances. Subsampling will occur once in every boosting iteration.
#' Range of (0, 1]. Defaults to 0.6.
#' @param min_child_weight minimum sum of instance weight required in a child node.
#' If the tree partitioning step produces a leaf node with a sum of instance weight
#'  less than \code{min_child_weight}, then the building process will cease further
#'   partitioning. A larger value of \code{min_child_weight} leads to a more
#'   conservative algorithm. Defaults to 1.
#' @param eta step size shrinkage. After each boosting step, one can obtain the
#' weights of new features directly, and the parameter \code{eta} is used to shrink
#' these feature weights, thereby making the boosting process more conservative.
#' Range of [0, 1]. Defaults to 0.3.
#' @param gamma minimum loss reduction needed to create an additional partition on
#' a leaf node of the tree. A larger value of \code{gamma} corresponds to a more
#' conservative algorithm. Defaults to 0.
#' @param max_delta_step maximum allowed step size for adjusting the output of each
#' leaf. If the value is set to 0, it indicates that there is no constraint. Defaults to 0.
#' @param lambda L2 regularization term on weights. Increasing this value will result in a more conservative model.
#' Defaults to 1.
#' @param alpha L1 regularization term on weights. Increasing this value will result in a more conservative model.
#' Defaults to 0.
#' @param na.rm if \code{TRUE}, observations with \code{NA} values are deleted
#' from the population and sample data. For the XGB procedure complete
#' observations are required. Defaults to \code{FALSE}.
#' @param ... additional parameters to be passed to \code{xgboost}.
#'
#' @return An object of class \code{xgb}, \code{emdi}, which encompasses point estimates,
#' uncertainty, and confidence intervals at the domain level, along with details regarding
#' the \code{xgb} model. Various generic functions such as \code{summary}, \code{estimators}
#' and \code{map_plot} are applicable to a model of the class \code{xgb}.
#' @references
#' Krennmair, P., & Schmid, T. (2022). Flexible Domain Prediction Using Mixed Effects
#' Random Forests. Journal of Royal Statistical Society: Series C (Applied Statistics),
#' Vol.71, No. 5, 1865–1894.\cr \cr
#' Merfeld, J. D., & Newhouse, D. (2023). Improving Estimates of Mean Welfare and Uncertainty
#' in Developing Countries (No. 10348). The World Bank.
#' @export
#' @importFrom xgboost xgboost
#' @importFrom dplyr select left_join arrange
#' @importFrom magrittr %>%
#' @importFrom purrr as_vector
#' @importFrom stats weighted.mean
#' @examples
#' \donttest{
#' # Loading data - population and sample data
#' data("eusilcA_pop")
#' data("eusilcA_smp")
#'
#' # Create subdomains: Equal to individuals in each area
#' eusilcA_smp$subDomain <- ave(eusilcA_smp$district,
#'                              eusilcA_smp$district,
#'                              FUN = seq_along)
#' eusilcA_pop$subDomain <- ave(eusilcA_pop$district,
#'                              eusilcA_pop$district,
#'                              FUN = seq_along)
#'
#' # Estimate extreme gradient boosting model
#' xgb_model <- xgb(fixed = eqIncome ~ eqsize + cash + self_empl +
#'                  unempl_ben + age_ben + surv_ben + sick_ben+
#'                  dis_ben +rent + fam_allow + house_allow +
#'                  cap_inv + tax_adj + district + subDomain,
#'                  smp_data = eusilcA_smp,
#'                  pop_data = eusilcA_pop,
#'                  domains = "district",
#'                  sub_domains = "subDomain")
#'
#' # Extract Mean, MSE and CV
#' estimators(object = xgb_model, indicator = "Mean",
#'            MSE = TRUE, CV =TRUE)
#'
#' # Plot the results on a map
#' load_shapeaustria()
#' map_plot(object = xgb_model, MSE = FALSE, CV = TRUE,
#'          map_obj = shape_austria_dis, indicator = c("Mean"),
#'          map_dom_id = "PB")
#'}

xgb <- function(fixed,
                smp_data,
                smp_weights = NULL,
                pop_data,
                pop_weights = NULL,
                domains,
                sub_domains,
                transformation = "no",
                B = 100,
                conf_level = 0.95,
                nround = 100,
                max_depth = 4,
                colsample_bytree = 0.6,
                colsample_bylevel = 0.6,
                colsample_bynode = 0.6,
                subsample = 0.6,
                min_child_weight = 1,
                eta = 0.3,
                gamma = 0,
                max_delta_step = 0,
                lambda = 1,
                alpha = 0,
                na.rm = FALSE,
                ...){

  out_call <- match.call()

  # Framework for xgb
  #_____________________________________________________________________________
  fwk <- framework_xgb(fixed = fixed,
                       smp_data = smp_data,
                       pop_data = pop_data,
                       smp_weights = smp_weights,
                       pop_weights = pop_weights,
                       domains = domains,
                       transformation = transformation,
                       conf_level = conf_level,
                       sub_domains = sub_domains,
                       na.rm = na.rm)

  # Direct estimates
  #_____________________________________________________________________________
  # Subdomains
  sub_domains_direct <- data.frame(cbind(fwk$Y_smp,
                                         fwk$X_smp[[paste0(sub_domains)]],
                                         fwk$X_smp[[paste0(domains)]]))
  colnames(sub_domains_direct) <- c("outcome", "sub_domains", "domains")
  sub_domains_direct$outcome <- as.numeric(sub_domains_direct$outcome)

  # Domains
  domains_direct <- data.frame(cbind(fwk$Y_smp,
                                     fwk$smp_weights,
                                     fwk$X_smp[[paste0(domains)]]))
  colnames(domains_direct) <- c("outcome", "wts", "domains")
  domains_direct$outcome <- as.numeric(domains_direct$outcome)
  domains_direct$wts <- as.numeric(domains_direct$wts)
  grouped_domains1 <- split(domains_direct$outcome, domains_direct$domains)
  weighted_means1 <- sapply(grouped_domains1, function(group) {
    stats::weighted.mean(group, wts = domains_direct$wts[domains_direct$domains == names(group)])
  })
  domains_direct <- data.frame(
    domains = names(weighted_means1),
    outcome = weighted_means1,
    row.names = NULL
  )

  # XGBoost
  #_____________________________________________________________________________
  X_smp_xgb <- fwk$X_smp %>%
    dplyr::select(-c(paste0(domains), paste0(sub_domains)))
  X_pop_xgb <- fwk$X_pop %>%
    dplyr::select(-c(paste0(domains), paste0(sub_domains)))

  xgb_fit <- xgboost::xgboost(
    data = as.matrix(X_smp_xgb),
    label = sub_domains_direct$outcome,
    weight = (fwk$smp_weights/mean(fwk$smp_weights)),
    nrounds = nround,
    max_depth = max_depth,
    colsample_bytree = colsample_bytree,
    colsample_bylevel = colsample_bylevel,
    colsample_bynode = colsample_bynode,
    subsample = subsample,
    min_child_weight = min_child_weight,
    eta = eta,
    gamma = gamma,
    max_delta_step = max_delta_step,
    lambda = lambda,
    alpha = alpha,
    verbose = 0,
    ...
  )

  # Predictions
  #_____________________________________________________________________________
  sub_pred <- data.frame(
    cbind(
      fwk$X_pop[[paste0(domains)]],
      fwk$X_pop[[paste0(sub_domains)]],
      fwk$pop_weights,
      predict(xgb_fit, as.matrix(X_pop_xgb))
    )
  )
  colnames(sub_pred) <- c("domains", "sub_domains", "wts", "hat")
  sub_pred$hat <- as.numeric(sub_pred$hat)
  sub_pred$wts <- as.numeric(sub_pred$wts)

  if (transformation=="arcsin"){
    sub_pred$hat <- ifelse(sub_pred$hat>asin(1), asin(1), sub_pred$hat)
    sub_pred$hat <- ifelse(sub_pred$hat<asin(0), asin(0), sub_pred$hat)
  }

  # Residuals
  #_____________________________________________________________________________
  # Subdomains
  sub_domains_direct <- sub_domains_direct %>%
    #dplyr::left_join(sub_pred, by = c("sub_domains", "domains"))
    dplyr::inner_join(sub_pred, by = c("sub_domains", "domains"))
  resid_sub_domains <- as.numeric(sub_domains_direct$outcome) - as.numeric(sub_domains_direct$hat)
  grouped_domains2 <- split(sub_pred$hat, sub_pred$domains)
  weighted_means2 <- sapply(grouped_domains2, function(group) {
    stats::weighted.mean(group, wts = as.numeric(sub_pred$wts[sub_pred$domains == names(group)]))
  })
  domains_pred <- data.frame(
    domains = names(weighted_means2),
    hat = weighted_means2,
    row.names = NULL
  )
  domains_direct <- domains_direct %>%
    dplyr::left_join(domains_pred, by = "domains")
  resid_domains <- as.numeric(domains_direct$outcome) - as.numeric(domains_direct$hat)

  # Bootstrap
  #_____________________________________________________________________________
  B_results <- matrix(data = NA, nrow = B, ncol = nrow(domains_pred))

  for (j in 1:B){

    B_sub <- sub_pred
    B_sub$hat <- as.numeric(B_sub$hat) + as.numeric(resid_sub_domains[sample(1:length(resid_sub_domains), nrow(B_sub), replace = TRUE)])
    grouped_domains3 <- split(B_sub$hat, B_sub$domains)
    weighted_means3 <- sapply(grouped_domains3, function(group) {
      stats::weighted.mean(as.numeric(group), wts = as.numeric(B_sub$wts[B_sub$domains == names(group)]))
    })
    B_domains <- data.frame(
      domains = names(weighted_means3),
      hat = weighted_means3,
      row.names = NULL
    )

    B_domains$hat <- as.numeric(B_domains$hat) + as.numeric(resid_domains[sample(1:length(resid_domains), nrow(B_domains), replace = TRUE)])
    B_domains <- B_domains %>%
      dplyr::arrange(domains)

    B_results[j,] <- purrr::as_vector(B_domains$hat)
  }

  # Prepare results
  #_____________________________________________________________________________
  sorted_results <- domains_pred[order(domains_pred$domains), ]
  results <- sorted_results[, c("domains", "hat")]

  results$lower <- NA
  results$upper <- NA
  results$sd <- NA

  for (l in 1:nrow(results)){

    temp <- B_results[,l]

    if (transformation=="arcsin"){
      temp <- ifelse(temp>asin(1), asin(1), temp)
      temp <- ifelse(temp<asin(0), asin(0), temp)
    }
    if (transformation=="arcsin"){
      temp <- sin(temp)^2
    }
    if (transformation=="log"){
      temp <- sin(exp(temp))^2
    }

    results$hat[l] <- mean(temp)
    results$lower[l] <- quantile(temp, probs = (1-conf_level)/2)
    results$upper[l] <- quantile(temp, probs = 1-(1-conf_level)/2)
    results$sd[l] <- sd(temp)
  }
  colnames(results) <- c("Domain", "Mean", "Lower", "Upper", "SD")


  result <- list(
    ind = data.frame(cbind(Domains = results["Domain"], Mean = results["Mean"])),
    MSE = data.frame(cbind(Domains = results["Domain"], Mean = results$SD)),
    CI  = data.frame(cbind(Domains = results["Domain"],
                           LowerCI = results["Lower"],
                           UpperCI = results["Upper"])),
    xgbModel = c(xgb_fit, call = out_call, smp_data = list(smp_data), transformation = transformation,
                 fwk$saeinfo)
  )
  class(result) <- c("xgb","emdi")
  return(result)

}




