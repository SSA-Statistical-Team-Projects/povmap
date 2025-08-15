#' Tuning extreme gradient boosting models for domain-level averages
#'
#' The funtion \code{xgb_tune} fine-tunes the hyperparameters for extreme gradient boosting models, following \cite{Merfeld and Newhouse (2023)}.
#' It offers the flexibility to allocate \code{domains} to folds and utilizes estimated
#' means at the domain-level for cross-validation. Users can specify the number of
#' folds, including an option for leave-one-out cross-validation.
#'
#' @param fixed a two-sided linear formula object describing the
#' fixed-effects part of the model with the dependent variable on the left
#' of a ~ operator and the explanatory variables on the right, separated
#' by + operators. All variables (except for \code{domains} and \code{cluster})
#' must be numeric.
#' @param smp_data a data frame that needs to comprise all variables including
#' \code{domains} and \code{cluster}.
#' @param smp_weights a character string containing the name of the variable that
#' indicates weights in \code{smp_data}. The variable has to be numeric.
#' Defaults to \code{NULL}.
#' @param domains a character string containing the name of a variable
#' that indicates domains in \code{smp_data}. The variable can be
#' numeric or a factor.
#' @param cluster a character string containing the name of a variable
#' that indicates clusters in \code{smp_data}. The variable can be
#' numeric or a factor. Defaults to \code{"domains"}.
#' @param transformation a character string. Two different transformation
#' types for the dependent variable can be chosen (i) no transformation ("no");
#' (ii) log transformation ("log"); (iii) Arcsin transformation ("arcsin").
#' Defaults to \code{"no"}.
#' @param folds number of folds. Defaults to 10.
#' @param nround combination of maximum number of boosting iterations. Defaults to 150 and 250.
#' @param max_depth combination of maximum depth of a tree. Defaults to 4 and 6.
#' @param colsample_bytree combination of subsample ratios of columns when constructing each tree. Defaults to 0.6 and 1.
#' @param colsample_bylevel combination of subsample ratios of columns for each level. Defaults to 0.6 and 1.
#' @param colsample_bynode combination of subsample ratios of columns for each node (split). Defaults to 0.6 and 1.
#' @param subsample combination of subsample ratios of the training instances. Defaults to 0.6 and 1.
#' @param min_child_weight minimum sum of instance weight required in a child node.
#' If the tree partitioning step produces a leaf node with a sum of instance weight
#' less than \code{min_child_weight}, then the building process will cease further
#' partitioning. A larger value of \code{min_child_weight} leads to a more conservative
#' algorithm. Defaults to 1.
#' @param eta step size shrinkage. After each boosting step, one can obtain the
#' weights of new features directly, and the parameter \code{eta} is used to shrink
#' these feature weights, thereby making the boosting process more conservative.
#' Range of [0, 1]. Defaults to 0.3.
#' @param gamma minimum loss reduction needed to create an additional partition
#' on a leaf node of the tree. A larger value of \code{gamma} corresponds to a more
#' conservative algorithm. Defaults to 0.
#' @param max_delta_step maximum allowed step size for adjusting the output of each
#' leaf. If the value is set to 0, it indicates that there is no constraint Defaults to 0.
#' @param lambda L2 regularization term on weights. Increasing this value will result in a more conservative model.
#' Defaults to 1.
#' @param alpha L1 regularization term on weights. Increasing this value will result in a more conservative model.
#' Defaults to 0.
#' @param verbose display progress. Defaults to FALSE.
#' @param ... additional parameters to be passed to \code{xgboost}.
#'
#' @return An object of class \code{xgb}, \code{emdi}, containing the optimal
#' hyperparameters for an extreme gradient boosting model.
#' @references
#' Merfeld, J. D., & Newhouse, D. (2023). Improving Estimates of Mean Welfare and Uncertainty
#' in Developing Countries (No. 10348). The World Bank. \cr \cr
#' @export
#' @importFrom xgboost xgboost
#' @importFrom dplyr left_join
#'
#' @examples
#' \donttest{
#' # Loading data - population and sample data
#' data("eusilcA_pop")
#' data("eusilcA_smp")
#'
#' xgb_tune_model <- xgb_tune(fixed = eqIncome ~ eqsize + cash + self_empl +
#'                            unempl_ben + age_ben + surv_ben + sick_ben +
#'                            dis_ben + rent + fam_allow + house_allow +
#'                            cap_inv + tax_adj + district,
#'                            smp_data = eusilcA_smp,
#'                            domains = "district")
#'}

xgb_tune <- function(fixed,
                     smp_data,
                     smp_weights = NULL,
                     domains,
                     cluster = "domains",
                     transformation = "no",
                     folds = 10,
                     nround = c(150, 250),
                     max_depth = c(4, 6),
                     colsample_bytree = c(0.6, 1),
                     colsample_bylevel = c(0.6, 1),
                     colsample_bynode = c(0.6, 1),
                     subsample = c(0.6, 1),
                     min_child_weight = c(1),
                     eta = c(0.3),
                     gamma = c(0),
                     max_delta_step = c(0),
                     lambda = c(1),
                     alpha = c(0),
                     verbose = FALSE,
                     ...){

  # Data preparation
  #_____________________________________________________________________________
  outcome <- all.vars(fixed[[2]])
  covariates <- all.vars(fixed[[3]])
  
  #split <- strsplit(as.character(fixed), "~", fixed = TRUE)
  #outcome <- trimws(split[[1]][1])
  #covariates <- trimws(strsplit(trimws(split[[1]][2]), "\\+")[[1]])
  X_smp <- smp_data[,c(covariates,domains)]
  Y_smp <- data.frame(smp_data[,outcome])

  if(is.null(smp_weights)==FALSE){
    smp_weights <- smp_data[,smp_weights]
  } else {
    smp_weights <- NULL
  }
  if (is.null(smp_weights)==TRUE){
    smp_weights <- rep(1, length = nrow(Y_smp))
  }
  if (cluster=="domains"){
    cluster <- paste0(domains)
  }
  colnames(Y_smp) <- "labels"

  # Check
  #_____________________________________________________________________________
  xgb_check2(
    transformation = transformation,
    Y_smp = Y_smp,
    X_smp = X_smp,
    smp_weights = smp_weights,
    domains = domains,
    cluster = cluster)

  # Transformation
  #_____________________________________________________________________________
  if (transformation=="arcsin"){
    Y_smp <- asin(sqrt(Y_smp))
  }
  if (transformation=="log"){
    Y_smp <- log(Y_smp)
  }

  # Folds
  #_____________________________________________________________________________
  cluster_col <- data.frame(X_smp[[paste0(cluster)]])
  colnames(cluster_col) <- paste0(cluster)
  cluster_unique <- data.frame(unique(cluster_col[,1]))
  colnames(cluster_unique) <- paste0(cluster)
  cluster_unique$fold <- sample(x = 1:folds, size = nrow(cluster_unique), replace = TRUE)

  cluster_col <- cluster_col %>%
    dplyr::left_join(cluster_unique, by = paste0(cluster))
  if (cluster=="domains"){
    X_final <- X_smp[,-c(which(colnames(X_smp)==paste0(cluster)))]
  } else{
    X_final <- X_smp[,-c(which(colnames(X_smp)==paste0(cluster)), which(colnames(X_smp)==paste0(domains)))]
  }

  # Grid
  #_____________________________________________________________________________
  tunegrid <- expand.grid(
    nround             = nround,
    max_depth          = max_depth,
    colsample_bytree   = colsample_bytree,
    colsample_bylevel  = colsample_bylevel,
    colsample_bynode   = colsample_bynode,
    subsample          = subsample,
    min_child_weight   = min_child_weight,
    eta                = eta,
    gamma              = gamma,
    max_delta_step     = max_delta_step,
    lambda             = lambda,
    alpha              = alpha
  )

  OPT <- matrix(NA, ncol = folds, nrow = dim(tunegrid)[1])

  # Tuning
  #_____________________________________________________________________________
  for (fold in 1:folds){

    for (row in 1:nrow(tunegrid)){

      xgb_fit <-  xgboost::xgboost(
        data               = data.matrix(X_final[cluster_col$fold!=fold,]),
        label              = Y_smp[cluster_col$fold!=fold,],
        weight             = smp_weights[cluster_col$fold!=fold],
        nrounds            = tunegrid$nround[row],
        max_depth          = tunegrid$max_depth[row],
        colsample_bytree   = tunegrid$colsample_bytree[row],
        colsample_bylevel  = tunegrid$colsample_bylevel[row],
        colsample_bynode   = tunegrid$colsample_bynode[row],
        subsample          = tunegrid$subsample[row],
        min_child_weight   = tunegrid$min_child_weight[row],
        eta                = tunegrid$eta[row],
        gamma              = tunegrid$gamma[row],
        max_delta_step     = tunegrid$max_delta_step[row],
        lambda             = tunegrid$lambda[row],
        alpha              = tunegrid$alpha[row],
        objective          = "reg:squarederror",
        verbose            = 0,
        ...
      )

      # Predictions (only for those out of sample)
      domains_hat <- data.frame(predict(xgb_fit, data.matrix(X_final[cluster_col$fold==fold,])))
      domains_hat[[paste0(domains)]] <- X_smp[cluster_col$fold==fold,][[paste0(domains)]]
      domains_hat[[colnames(Y_smp)]] <- Y_smp[cluster_col$fold==fold,]
      colnames(domains_hat) <- c("hat", "domains", "labels")
      grouped_domains <- split(domains_hat, domains_hat$domains)
      mean_hat <- sapply(grouped_domains, function(group) mean(group$hat))
      mean_labels <- sapply(grouped_domains, function(group) mean(group$labels))
      first_rows <- lapply(grouped_domains, function(group) group[1, ])
      domains_pred <- do.call(rbind, first_rows)
      domains_pred$hat <- mean_hat
      domains_pred$labels <- mean_labels

      # Predict, square error, and take mean --> MSE
      OPT[row, fold] <- mean((domains_pred$labels - domains_pred$hat)^2)
      if (verbose==TRUE){
        print(paste0("Fold ", fold, " of ", folds, " and row ", row, " of ", nrow(tunegrid)))
      }
    }
  }

  # Optimal values
  #_____________________________________________________________________________
  mse_min <- min(apply(OPT, 1, FUN = mean))
  nround_opt <- tunegrid$nround[which.min(mse_min)]
  max_depth_opt <- tunegrid$max_depth[which.min(mse_min)]
  colsample_bytree_opt <- tunegrid$colsample_bytree[which.min(mse_min)]
  colsample_bylevel_opt <- tunegrid$colsample_bylevel[which.min(mse_min)]
  colsample_bynode_opt <- tunegrid$colsample_bynode[which.min(mse_min)]
  subsample_opt <- tunegrid$subsample[which.min(mse_min)]
  min_child_weight_opt <- tunegrid$min_child_weight[which.min(mse_min)]
  eta_opt <- tunegrid$eta[which.min(mse_min)]
  max_delta_step_opt <- tunegrid$max_delta_step[which.min(mse_min)]
  gamma_opt <- tunegrid$gamma[which.min(mse_min)]
  lambda_opt <- tunegrid$lambda[which.min(mse_min)]
  alpha_opt <- tunegrid$alpha[which.min(mse_min)]

  final_output <- list(nround_opt, max_depth_opt, colsample_bytree_opt, colsample_bylevel_opt,
                       colsample_bynode_opt, subsample_opt, min_child_weight_opt, eta_opt,
                       gamma_opt, max_delta_step_opt, lambda_opt, alpha_opt)
  names(final_output) <- c("nround", "max_depth", "colsample_bytree", "colsample_bylevel",
                           'colsample_bynode', "subsample", "min_child_weight", "eta",
                           "gamma", "max_delta_step", 'lambda', "alpha")

  class(final_output) <- c("xgb","emdi")
  return(final_output)
}


