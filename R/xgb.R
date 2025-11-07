#' Extreme gradient boosting for domain-level averages
#'
#' The function \code{xgb} employs extreme gradient boosting to estimate domain-level averages, particularly
#' for small area estimation (SAE) applications. The model is specified 
#' at the sub-area level (any geopgraphic level more disaggregated than the target areas), 
#' as implemented by \cite{Merfeld, Dang, and Newhouse (2025)} To estimate the mean squared
#' error (MSE), a nonparametric residual bootstrap approach is utilized, as described
#' in \cite{Krennmair and Schmid (2022)} and \cite{Merfeld, Dang, and Newhouse (2025)}.
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
#' @param bootstrap If TRUE, implements bootstrap procedure to estimate variance. 
#' Defaults to TRUE. 
#' @param B a number determining the number of bootstrap populations in the
#' nonparametric residual bootstrap approach used in the MSE estimation. The
#' number must be greater than 1. Defaults to 1000. For practical applications,
#' values larger than 200 are recommended.
#' @param conf_level confidence level for the confidence interval. Defaults to 0.95.
#' @param nrounds maximum number of boosting iterations. Defaults to 100.
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
#' @param seed an integer to set the seed for the random number generator, see Details.
#' #' Defaults to 123
#' @param ydump if not null, specifies a file name to save XGboost predictions. Defaults to NULL. 
#' @param benchmark The input depends on the type of benchmarking to be 
#' performed.
#' (i) Benchmarking with a fixed value:
#' (a) with one value for each indicator: a named vector containing the numeric
#' benchmark value(s). The names of the vector matchs to the chosen indicators.
#' Benchmarking is available for \code{"Mean"} and \code{"Head_Count"}.
#' (b) with values for the sub-level specified in the argument
#' @param benchmark_level: a data.frame composed of a variable of class
#' character containing the domain names at which the benchmarking is
#' performed and variable(s) with benchmark value(s) of class numeric.
#' Benchmarking is supplied for the Mean and the Head_Count ratio. Therefore,
#' the names of the data.frame must match for the first variable the
#' benchmark_level and for the other(s) to Mean and Head_Count.
#' (ii) Benchmarking with the survey data: a vector containing the names of the
#' chosen indicators. In this case, survey weights (\code{weights}) are needed.
#' Benchmarking is available for \code{"Mean"} and \code{"Head_Count"}.
#' @param benchmark_type a character indicating the type of benchmarking. Types
#' that can be chosen (i) Raking ("\code{raking}"), (ii) Ratio adjustment
#' ("\code{ratio}"), (iii) ratio adjustment of the complement
#' ("\code{ratio_complement}" and (iv) ratio adjustment when the maximum 
#' benchmarked estimate in a benchmark level <=1 and ratio adjustment of the 
#' complement when the maximum benchmarked estimate >1 ("\code{ratio_bound}. 
#' Defaults to "\code{ratio}"
#' @param benchmark_level a character indicating the level at which the
#' benchmarking is performed. This name must be represented in the sample and
#' population data as a variable name.
#' @param benchmark_weights the name of variable containing benchmark weights.
#' This is only possible for internal benchmarking and enable users to benchmark
#' with weights differing from the survey weights (Default for weighting for
#' internal benchmarking).
#' @param weightedBS. If TRUE and smp_weights is specified, gives each area and subarea weight proportional to 
#' their sample weight when implementing the cluster residual bootstrap. If FALSE, gives each area and subarea 
#' equal weight. Defaults to TRUE.
#' @param center_residuals. If TRUE, residuals are centered around zero prior to adding then to XGboost predictions  
#' to generate estimates. Defaults to FALSE.    
#' @param ... additional parameters to be passed to \code{xgboost}.
#'
#' @return An object of class \code{xgb}, \code{emdi}, which includes point estimates,
#' uncertainty, and confidence intervals at the domain level, along with details regarding
#' the \code{xgb} model. Various generic functions such as \code{summary}, \code{estimators}
#' and \code{map_plot} are applicable to a model of the class \code{xgb}.
#' @references
#' Krennmair, P., & Schmid, T. (2022). Flexible Domain Prediction Using Mixed Effects
#' Random Forests. Journal of Royal Statistical Society: Series C (Applied Statistics),
#' Vol.71, No. 5, 1865–1894.\cr \cr
#' Merfeld, J. D., Dang, H., & Newhouse, D. (2025). Improving Estimates of Mean Welfare and Uncertainty
#' in Developing Countries (No. 10348). The World Bank.
#' @importFrom purrr as_vector
#' @importFrom collapse fmean 
#' @export
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
                bootstrap = T, 
                B = 1000,
                conf_level = 0.95,
                # use MORE CONSERVATIVE XGBoost defaults (https://xgboost.readthedocs.io/en/stable/parameter.html)
                nrounds = 100,
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
                seed = 123,
                benchmark = NULL,
                benchmark_type = "ratio",
                benchmark_level = NULL,
                benchmark_weights = NULL,
                ydump = NULL,
                weightedBS = T, 
                center_residuals = F, 
                ...){

  out_call <- match.call()
  
  if (is.null(benchmark_weights) & !is.null(smp_weights)) {
    benchmark_weights <- smp_weights
  }
  
  
  
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
                       na.rm = na.rm,
                       benchmark = benchmark, 
                       benchmark_level = benchmark_level,
                       benchmark_weights = benchmark_weights,
                       benchmark_type = benchmark_type)

  # Direct estimates
  #_____________________________________________________________________________
  # Subdomains
  sub_domains_direct <- data.frame(fwk$Y_smp,
                                         fwk$X_smp[sub_domains],
                                         fwk$X_smp[domains],
                                         fwk$smp_weights_vec)
  colnames(sub_domains_direct) <- c("outcome", "sub_domains", "domains","smp_weights")
  #sub_domains_direct$domains <- as.character(sub_domains_direct$domains)
  #sub_domains_direct$sub_domains <- as.character(sub_domains_direct$sub_domains)
  

  
  # Estimate XGBoost
  

  X_smp_xgb <- fwk$X_smp[c(fwk$covariates)]
  #X_smp_xgb[,c(paste0(domains),paste0(sub_domains))] <- list(NULL)
  
  #We only want the predictors 
  X_pop_xgb <- fwk$X_pop[,fwk$covariates]
  
  set.seed(seed)
  xgb_fit <- xgboost::xgboost(
    data = as.matrix(X_smp_xgb),
    label = sub_domains_direct$outcome,
    weight = (fwk$smp_weights_vec/mean(fwk$smp_weights_vec)),
    nrounds = nrounds,
    max_depth = max_depth, 
    colsample_bytree = colsample_bytree, colsample_bylevel = colsample_bylevel, 
    colsample_bynode = colsample_bynode, subsample = subsample, 
    min_child_weight = min_child_weight, eta = eta, gamma = gamma, max_delta_step = max_delta_step,
    lambda = lambda, alpha = alpha,
    verbose = 0, 
    nthread=1,
    ...
  )
  
  
  #Generate predictions at sub-area level (assumed to be each observation)
  
  sub_pred <- data.frame(
      fwk$X_pop[,domains],
      fwk$X_pop[,sub_domains],
      fwk$X_pop[,pop_weights],
      predict(xgb_fit, as.matrix(X_pop_xgb))
    )

  colnames(sub_pred) <- c("domains", "sub_domains", "wts", "hat")

  if (transformation=="arcsin"){
    sub_pred$hat <- ifelse(sub_pred$hat>asin(1), asin(1), sub_pred$hat)
    sub_pred$hat <- ifelse(sub_pred$hat<asin(0), asin(0), sub_pred$hat)
  }

  ### write predictions to file if needed
  if (!is.null(ydump)) {
    saveRDS(sub_pred, ydump)
  }

  
  
  
  
  
  

  # Residuals
  #_____________________________________________________________________________
  # Subdomains
  
  sub_domains_direct <-merge(x = sub_domains_direct, y = sub_pred, by = c("sub_domains", "domains"), all.x = T)
  sub_domains_direct <- sub_domains_direct[order(sub_domains_direct$sub_domains),]
  
  
  

  #sub_domains_direct <- sub_domains_direct[order(sub_domains_direct$sub_domains),]

  domains_pred <- aggregate_weighted_mean(df=sub_pred$hat,by=list(sub_pred$domains),w=sub_pred$wts) 
  # add sum of weights to domain-level predictions df
  wts <- aggregate(x=sub_pred$wts,by=list(sub_pred$domains),FUN = sum)
  wts[,1] <- unique(sub_pred$domains)
  domains_pred <- data.frame(domains_pred,wts$x)
  colnames(domains_pred) <- c("domains","hat","weight")
  
  # Get direct estimates at domain level 
  domains_direct <- data.frame(fwk$Y_smp,
                               fwk$smp_weights_vec,
                               fwk$X_smp[paste0(domains)])
  colnames(domains_direct) <- c("outcome", "wts", "domains")
  domains_direct <- aggregate_weighted_mean(df=domains_direct$outcome,by=list(domains_direct$domains),w=domains_direct$wts)
  colnames(domains_direct) <- c("domains","outcome")
  
  # benchmarking exercise to obtain modified residuals for bootstrap (benchmarking of points estimates occurs at end of function below)
  #domains_pred$hat_bench <- domains_pred$hat 
  #sub_domains_direct$hat_bench <- sub_domains_direct$hat 
  # if (!is.null(benchmark)) {
  #   point_estim <- NULL
  #   point_estim$ind <- data.frame(Mean = domains_pred$hat)
  #   if (is.null(benchmark_level)) {
  #     point_estim$ind <- benchmark_ebp_national(
  #       point_estim = point_estim,
  #       framework = fwk,
  #       fixed = fixed,
  #       benchmark = benchmark,
  #       benchmark_type = benchmark_type)
  #   } else {
  #     point_estim$ind <- benchmark_ebp_level(
  #       point_estim = point_estim,
  #       framework = fwk,
  #       fixed = fixed,
  #       benchmark = benchmark,
  #       benchmark_type = benchmark_type,
  #       benchmark_level = benchmark_level)
  #   }
  #   domains_pred$hat_bench <- point_estim$ind$Mean_bench
  #   }
  # adjust sample sub_domain predictions according to benchmark ratios  
  #benchmark_factors <- data.frame(domains = domains_pred$domains,benchmark_factor = domains_pred$hat_bench/domains_pred$hat)
  # benchmark_factor is a vector of ones when benchmarking not selected 
  #sub_domains_direct <- merge(sub_domains_direct,benchmark_factors,by="domains",all.X=T )
  resid_sub_domains <- (sub_domains_direct$outcome - as.numeric(sub_domains_direct$hat)) 
  
  
  # merge direct estimates at area level with predictions while maintaining sort order
  domains_direct <- merge(x=domains_direct,y=domains_pred,by="domains",all.x=T)
  #domains_direct <- merge(x=domains_direct,y=benchmark_factors,by="domains",all.x=T)

  # calculate domain residuals 
    resid_domains <- (domains_direct$outcome - domains_direct$hat) 


  
  
  #pop_subarea_d <- (subarea_coords$weights*subarea_coords$n)/(sum(subarea_coords$weights*subarea_coords$n))
  #pop_area_d <- (area_coords$weights*area_coords$n)/(sum(area_coords$weights*area_coords$n))

  # Bootstrap
  #_____________________________________________________________________________
  B_results <- matrix(data = NA, nrow = B, ncol = nrow(domains_pred))
  B_results_bench <- matrix(data = NA, nrow = B, ncol = nrow(domains_pred))
  colnames(B_results) <- unique(sub_pred$domains)
  
  # Sample with probability proportional to sample weight if weightsBS = TRUE
  if (weightedBS==T & !is.null(smp_weights)) {
    pop_subarea_d <- sub_domains_direct$smp_weights/sum(sub_domains_direct$smp_weights)  
    domain_wts <- aggregate(sub_domains_direct$smp_weights,by=list(sub_domains_direct$domains),FUN=sum)$x
    pop_area_d <- domain_wts/sum(domain_wts)
  } else {
    # Otherwise sample each subarea and area with fixed proobability 
    pop_subarea_d <- rep(1/length(resid_sub_domains),length(resid_sub_domains))
    pop_area_d <- rep(1/length(resid_domains),length(resid_domains))
  }
  B_sub <- sub_pred  
  B_sub <- B_sub[order(B_sub$sub_domains), ]
  results <- domains_pred[, c("domains", "hat")]
  
  
if (bootstrap==T) {  
cat("Beginning bootstrap \n")
  for (j in 1:B){
    displayevery = max(round(B/10,0),1)
    
    if (j %% displayevery==0) {
      cat(paste0("replication ",j," of ",B,"\n"))
    }

    #bs_sub <- bs_sub |>
    #  arrange(subarea)
    # randomly sample residuals USING THE WEIGHTS CALCULATED ABOVE if weighted_BS==TRUE
    B_sub$sim <- as.numeric(B_sub$hat) + resid_sub_domains[sample(1:length(resid_sub_domains),
                                                                             nrow(B_sub), prob=pop_subarea_d, 
                                                                             replace = TRUE)]
    
    # aggregate to area
    #B_domains <- aggregate_weighted_mean(df=B_sub$sim,by=list(B_sub$domains),w=B_sub$wts)
    #colnames(B_domains) <- c("domains","sim")  
    B_domains <- collapse:::fmean(x=B_sub$sim,g=B_sub$domains,w=B_sub$wts)
    B_domains <- data.frame("domains" = names(B_domains),"sim" = B_domains)
                 
    area_draws <- data.frame(domains=B_domains$domains,area_draw=resid_domains[sample(1:length(resid_domains),
                                                                                      nrow(B_domains), prob=pop_area_d, 
                                                                                      replace = TRUE)]
    )
    # randomly sample residuals USING THE WEIGHTS CALCULATED ABOVE if weighted_BS==TRUE
    B_domains$sim <- B_domains$sim +area_draws$area_draw     
    
    
    if (!is.null(benchmark)) {
      B_sub[,sub_domains] <- B_sub$sub_domains
    B_sample <- merge(smp_data,B_sub,by=sub_domains,all.x=T)[,c(sub_domains,"sim","domains",benchmark_level,benchmark_weights)]
    B_sample_d <- collapse:::fmean(B_sample$sim,g=B_sample$domains,w=B_sample[,benchmark_weights])
    B_sample_d <- data.frame("domains" = names(B_sample_d),"sim" = B_sample_d)
    B_sample_d <- data.frame(B_sample_d, collapse:::ffirst(B_sample[,benchmark_level],g=B_sample$domains))
    B_sample_d <- data.frame(B_sample_d, collapse:::fsum(B_sample[,benchmark_weights],g=B_sample$domains))
    colnames(B_sample_d)[3:4] <- c(benchmark_level,benchmark_weights) 
    
    # create new area draws 
    B_sample_d <- merge(B_sample_d,area_draws,by="domains",all.X=T)
    B_sample_d$sim <- B_sample_d$sim+B_sample_d$area_draw
    
#    area_draws2 <- data.frame(domains=B_sample_d$domains,area_draw=resid_domains[sample(1:length(resid_domains),
#                                                                                      nrow(B_sample_d), prob=pop_area_d, 
#                                                                                      replace = TRUE)])
 #   B_sample_d$sim <- B_sample_d$sim+area_draws2$area_draw
    
    bm <- collapse:::fmean(B_sample_d$sim,g=B_sample_d[,benchmark_level],w=B_sample_d[,benchmark_weights])
    bm <- data.frame("state_" = names(bm),"Mean" = bm)
    point_estim <- NULL
    point_estim$ind <- data.frame("Mean" = B_domains$sim)
    if (is.null(benchmark_level)) {
          point_estim$ind <- benchmark_ebp_national(
            point_estim = point_estim,
            framework = fwk,
            fixed = fixed,
            benchmark = "Mean",
            benchmark_type = benchmark_type)
        } else {
          point_estim$ind <- benchmark_xgb_level(
            point_estim = point_estim,
            framework = fwk,
            fixed = fixed,
            benchmark = bm,
            benchmark_type = benchmark_type,
            benchmark_level = "state_")
        }
    B_results_bench[j,] <- purrr::as_vector(point_estim$ind$Mean_bench)
    } # close benchmarking loop for bootstrap 
    
    
    
    B_results[j,] <- purrr::as_vector(B_domains$sim)
} # close bootstrap loop 
  # Prepare results
  #_____________________________________________________________________________

  results$lower <- NA
  results$upper <- NA
  results$var <- NA
  results$var_bench <- NULL 
  
  for (l in 1:nrow(results)){

    temp <- B_results[,l]
   
    if (center_residuals==T) {
            temp <- temp + domains_pred$hat[l] - mean(temp) 
    }
    
    if (transformation=="arcsin"){
      temp <- ifelse(temp>asin(1), asin(1), temp)
      temp <- ifelse(temp<asin(0), asin(0), temp)
    }
    if (transformation=="arcsin"){
      temp <- sin(temp)^2
    }
    if (transformation=="log"){
      temp <- exp(temp)
    }
        
    results$hat[l] <- mean(temp)
    results$lower[l] <- quantile(temp, probs = (1-conf_level)/2)
    results$upper[l] <- quantile(temp, probs = 1-(1-conf_level)/2)
    results$var[l] <- var(temp)
    if (!is.null(benchmark)) {
      temp_bench <- B_results_bench[,l]
      results$var_bench[l] <- var(temp_bench)
    }
    
  }

  colnames(results) <- c("Domain", "Mean", "Lower", "Upper", "var")
  if (!is.null(benchmark)) {
  colnames(results)[6] <- c("var_bench")
    }
  

  result <- list(
    ind = data.frame(Domain = results["Domain"], Mean = results["Mean"]),
    var = data.frame(Domain = results["Domain"], Mean = results["var"]),
    CI  = data.frame(Domain = results["Domain"],
                           LowerCI = results["Lower"],
                           UpperCI = results["Upper"]),
    yhat=data.frame(sub_domains_direct[,c("sub_domains","hat")]),
    model = xgb_fit, 
    smp_data =  smp_data, 
    out_call = out_call, 
    transformation = transformation, 
    framework = fwk
  )
  colnames(result$var)[2]="Mean"
  } else {
    # Bootstrap not selected 
    result <- list(
    ind = data.frame(Domain = domains_pred[,"domains"], Mean = domains_pred[,"hat"]),
    yhat=data.frame(sub_domains_direct[,c("sub_domains","hat")]),
    model = xgb_fit, 
    smp_data =  smp_data, 
    out_call = out_call, 
    transformation = transformation, 
    framework = fwk
    ) 
  }

  if (!is.null(benchmark)) {
    point_estim <- NULL
    point_estim$ind <- data.frame(Mean = result$ind$Mean)
    
    if (is.null(benchmark_level)) {
      point_estim$ind <- benchmark_ebp_national(
        point_estim = point_estim,
        framework = fwk,
        fixed = fixed,
        benchmark = benchmark,
        benchmark_type = benchmark_type)
    } else {
      point_estim$ind <- benchmark_xgb_level(
        point_estim = point_estim,
        framework = fwk,
        fixed = fixed,
        benchmark = benchmark,
        benchmark_type = benchmark_type,
        benchmark_level = benchmark_level)
    }
    result$ind$Mean_bench <- point_estim$ind$Mean_bench
    colnames(result$ind)[ncol(result$ind)] <- "Mean_bench"
    result$var$Mean_bench <- results$var_bench
    colnames(result$var)[ncol(result$var)] <- "Mean_bench"
  }
  
  class(result) <- c("xgb","povmap")
  return(result)
}
