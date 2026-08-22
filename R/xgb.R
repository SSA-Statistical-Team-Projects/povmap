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
#' @param transformation a character string. Transformation
#' types for the dependent variable: (i) no transformation ("no");
#' (ii) log transformation ("log"); (iii) Arcsin transformation ("arcsin");
#' (iv) log-shift transformation ("log.shift"); (v) logistic transformation ("logistic");
#' (vi) Poisson log1p transformation for count data ("poisson").
#' Defaults to \code{"no"}.
#' @param bootstrap If TRUE, implements bootstrap procedure to estimate variance.
#' Defaults to TRUE.
#' @param bootstrap_type character spring specifying the type of bootstrap implemented.
#' (i) "case" requests a cluster case resampling bootstrap while (ii) "residual" requests a
#' cluster residual bootstrap. Defaults to "residual".
#' @param B a number determining the number of bootstrap populations in the
#' nonparametric residual bootstrap approach used in the MSE estimation. The
#' number must be greater than 1. Defaults to 200. For practical applications,
#' values larger than 200 are recommended.
#' @param L a number determining the number of Monte-Carlo residual bootstrap simulations to run to generate
#' point estimates when using a transformation. Defaults to 100.
#' @param weightedBS. If TRUE and smp_weights is specified, gives each area and subarea weight proportional to
#' their sample weight when implementing the cluster residual bootstrap. If FALSE, gives each area and subarea
#' equal weight. Defaults to TRUE.
#' @param boot_estimates. If TRUE, point_estimates are set equal to the average of the bootstrap replications. If set to FALSE, point estimates
#' are set equal to the XGboost prediction. Defaults to FALSE.
#' @param center_residuals. If TRUE, the mean of the residuals during the bootstrap are subtracted from the residuals prior to the bootstrapping procedure.
#' Defaults to FALSE.
#' @param cpus. Number of cores to parallelize across. Defaults to 1 (no parallelization)
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
#' @param variance_y the name of the variable in the sample data that contains the variance of the outcome variable.
#' This is used to account for heteroscedasticity when estimating the model. If NULL, no adjustment is made for heteroscedasticity.
#' Defaults to NULL.
#' @param rescale_weights If TRUE, rescales sample weights within each domain so that
#' they sum to the domain sample size. This prevents populous domains from dominating
#' the loss function, analogous to the treatment of weights in mixed models.
#' Defaults to TRUE.
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
#' @param perturb_benchmark logical. If \code{TRUE} and \code{benchmark_level}
#'   is specified, the bootstrap benchmark target at each benchmark-level group
#'   (e.g. state) is perturbed in every bootstrap iteration by adding
#'   independent Gaussian noise whose standard deviation equals the
#'   Horvitz-Thompson estimate of the direct-estimate standard error at that
#'   group. This propagates the survey-sampling uncertainty of the benchmark
#'   target into the benchmarked confidence intervals, correcting potential
#'   undercoverage that arises when the bootstrap state mean is far more
#'   stable than the real direct estimate. Defaults to \code{FALSE}.
#' @param benchmark_target_se optional named numeric vector giving the standard
#'   error of each benchmark-level direct estimate, keyed by the
#'   \code{benchmark_level} group value. When supplied (and
#'   \code{perturb_benchmark = TRUE}), it is used as the SD of the
#'   benchmark-target perturbation instead of the internal Horvitz-Thompson
#'   estimate. Supply the same estimator used to build the published direct-
#'   estimate confidence intervals (typically the cluster-robust / Taylor SE) so
#'   the perturbation injects exactly the uncertainty the target is reported
#'   with. An error is raised if any benchmark group is missing from the vector.
#'   Defaults to \code{NULL} (use the internal Horvitz-Thompson SE).
#' @importFrom purrr as_vector
#' @importFrom collapse fmean
#' @importFrom foreach foreach %dopar% %do%
#' @importFrom doParallel registerDoParallel
#' @importFrom doSNOW registerDoSNOW txtProgressBar
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

# Horvitz-Thompson variance of a weighted mean for each group.
# Returns a named numeric vector of variances (one entry per group).
# Formula: V_HT(ȳ_g) = n_g / ((n_g-1) * (Σw)²) * Σ w²(y - ȳ)²
# Equivalent to the with-replacement linearisation used by survey::svymean.
ht_var_weighted_mean <- function(y, w, g) {
  # Horvitz-Thompson / Hajek variance of a weighted MEAN (ratio estimator
  # yhat = sum(w y) / sum(w)) under Poisson sampling:
  #   sigma_hat^2 = (1 / (sum w)^2) * sum_i w_i (w_i - 1) (y_i - yhat)^2
  # The residual (y_i - yhat) is REQUIRED: sum w(w-1) y_i^2 is the variance of
  # the TOTAL (sum w y); for the ratio mean the deviations from the estimated
  # mean must be used, otherwise the variance scales with E[y^2] rather than the
  # dispersion of y and is grossly inflated for outcomes whose mean is far from
  # zero (e.g. a proportion near 0.6: ~7x too large in variance, ~3x in SE).
  # Census-validated: the centered form reproduces the design-based province
  # cluster-robust (Taylor) SE used elsewhere in the project (0.043 vs 0.044 for
  # DRC poverty); the uncentered form gave 0.144.
  g <- as.character(g)
  groups <- unique(g)
  result <- setNames(numeric(length(groups)), groups)
  for (grp in groups) {
    idx  <- g == grp
    y_g  <- y[idx];  w_g <- w[idx];  n_g <- sum(idx)
    if (n_g < 2L) { result[grp] <- NA_real_; next }
    ybar_g <- sum(w_g * y_g) / sum(w_g)
    result[grp] <- sum(w_g * (w_g - 1) * (y_g - ybar_g)^2) / sum(w_g)^2
  }
  result
}

# Tested-version notice.
# The xgb point estimate and bootstrap are sensitive to the xgboost version:
# cross-version prediction drift is diffuse (~0.05 median per ward between
# 1.7.7.1 and 3.1.2.1 on identical data) and is NOT fixable by setting
# base_score/tree_method/max_bin explicitly. povmap's results were validated
# against xgboost 3.1.2.1, so estimates produced on another version may differ
# from published ones even on identical data and tuning.
#
# This used to stop(). It no longer does: pinning to a single version makes the
# package unusable on every future xgboost release, which is a worse problem than
# the drift it guards against. The finding is preserved as a one-time warning so
# the user can judge whether cross-version comparability matters for their use.
# Suppress with options(povmap.skip_xgb_version_check = TRUE).
.povmap_xgb_tested <- "3.1.2.1"
.povmap_xgb_warned <- new.env(parent = emptyenv())
.assert_xgb_version <- function() {
  if (isTRUE(getOption("povmap.skip_xgb_version_check", FALSE))) return(invisible(NULL))
  found <- as.character(utils::packageVersion("xgboost"))
  if (identical(found, .povmap_xgb_tested)) return(invisible(NULL))
  # warn once per session: three entry points call this, and xgb_tune/xgb_cv call
  # it inside loops, so a per-call warning would bury the message it is making.
  if (isTRUE(.povmap_xgb_warned$done)) return(invisible(NULL))
  .povmap_xgb_warned$done <- TRUE
  warning(sprintf(
    paste0("xgboost %s detected; povmap's xgb results were validated against %s.\n",
           "  Point estimates and bootstrap intervals are version-sensitive: drift of\n",
           "  ~0.05 median per ward was observed between 1.7.7.1 and 3.1.2.1 on identical\n",
           "  data and tuning, and is not removable by setting base_score, tree_method or\n",
           "  max_bin. Results remain internally consistent on any single version; only\n",
           "  cross-version comparisons are affected.\n",
           "  (Silence with options(povmap.skip_xgb_version_check = TRUE).)"),
    found, .povmap_xgb_tested), call. = FALSE)

  invisible(NULL)
}

xgb <- function(fixed,
                smp_data,
                smp_weights = NULL,
                pop_data,
                pop_weights = NULL,
                domains,
                sub_domains,
                transformation = "no",
                bootstrap = T,
                bootstrap_type = "residual",
                B = 200,
                L=100,
                cpus=1,
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
                boot_estimates = F,
                center_residuals = F,
                variance_y = NULL,
                rescale_weights = TRUE,
                perturb_benchmark = FALSE,
                benchmark_target_se = NULL,
                verbose = FALSE,
                ...){

  #1. Initialize
  .assert_xgb_version()
  out_call <- match.call()
  # default to using sample weights for benchmarking if internal benchmarking
  if (is.null(benchmark_weights) & !is.null(smp_weights)) {
    benchmark_weights <- smp_weights
  }
  collapse:::set_collapse(sort = FALSE)

  # 1. Framework for xgb
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
                       benchmark_type = benchmark_type,
                       variance_y = variance_y)


  # 2. Obtain direct estimates
  #_____________________________________________________________________________
  #Transform outcome if called for
  if (transformation=="arcsin"){
    transform_outcome <- arcsin_transform
    back_transform_outcome <- arcsin_transform_back
  }
  else if (transformation=="log"){
    transform_outcome <- log_transform
    back_transform_outcome <- log_transform_back
  }
  else if (transformation=="log.shift") {
    # Compute shift from sample data: ensures min(y) + ls_lambda > 0
    ls_lambda <- if (min(fwk$Y_smp) <= 0) abs(min(fwk$Y_smp)) + 1 else 0
    transform_outcome      <- function(y) list(y = log(y + ls_lambda), shift = NULL)
    back_transform_outcome <- function(y) exp(y) - ls_lambda
  }
  else if (transformation=="logistic") {
    transform_outcome <- logit_transform_epsilon
    back_transform_outcome <- logit_transform_back
  }

  else if (transformation=="poisson") {
    if (min(fwk$Y_smp) < 0) stop("Outcome must be non-negative when using the poisson transformation")
    transform_outcome <- poisson_transform
    back_transform_outcome <- poisson_transform_back
  }
  else if (transformation=="no") {
    transform_outcome <- no_transform
    back_transform_outcome <- no_transform_back
    # Have to distinguish between no_transform and from no_transform_back because of shift parameter
  }
  else {
    stop("transformation must be 'no', 'log', 'log.shift', 'logistic', 'arcsin', or 'poisson'")
  }

  # Create a dataframe with sub_domains data
  sub_domains_direct <- data.frame(fwk$Y_smp,
                                   fwk$X_smp[sub_domains],
                                   fwk$X_smp[domains],
                                   fwk$smp_weights_vec)
  colnames(sub_domains_direct) <- c("outcome", sub_domains, "domains","smp_weights")


  # 3. Set up data to estimate XGBoost and generate point estimates

  X_smp_xgb <- fwk$X_smp[c(fwk$covariates)]
  #X_smp_xgb[,c(paste0(domains),paste0(sub_domains))] <- list(NULL)

  #We only want the predictors
  X_pop_xgb <- fwk$X_pop[,fwk$covariates]

  set.seed(seed)
  params <- list(max_depth          = max_depth,
                 colsample_bytree   = colsample_bytree,
                 colsample_bylevel  = colsample_bylevel,
                 colsample_bynode   = colsample_bynode,
                 subsample          = subsample,
                 min_child_weight   = min_child_weight,
                 eta                = eta,
                 gamma              = gamma,
                 max_delta_step     = max_delta_step,
                 lambda             = lambda,
                 alpha              = alpha,
                 nthread            = 1,
                 ...)

  # Rescale weights within each domain so they sum to the domain sample size
  # This prevents populous domains from dominating the loss function
  if (rescale_weights) {
    domain_sum_wts <- ave(fwk$smp_weights_vec, sub_domains_direct$domains, FUN = sum)
    domain_n <- ave(fwk$smp_weights_vec, sub_domains_direct$domains, FUN = length)
    smp_weights_rescaled <- fwk$smp_weights_vec * domain_n / domain_sum_wts
  } else {
    smp_weights_rescaled <- fwk$smp_weights_vec
  }

  # Estimate model
  xgb_model <- point_estim_xgb(params=params,smp_X=X_smp_xgb,
                               smp_Y=sub_domains_direct$outcome,
                               smp_weight= smp_weights_rescaled/mean(smp_weights_rescaled),
                               pop_X=X_pop_xgb,
                               sub_domains=sub_domains,
                               nrounds=nrounds,
                               transform_outcome=transform_outcome,
                               back_transform_outcome = back_transform_outcome,
                               L=L,
                               fwk=fwk,
                               variance_y = variance_y)

  xgb_fit <-  xgb_model$model
  domains_pred <- xgb_model$predictions

  ### write predictions to file if needed
  if (!is.null(ydump)) {
    saveRDS(xgb_model$sub_predictions, ydump)
  }

  resid_sub_domains <- xgb_model$resid_sub_domains
  resid_domains <- xgb_model$resid_domains



  #browser()
  #4. benchmark area-level estimates and transform these if option is selected, so bootstrap can draw from the residuals of the benchmarked estimates
  if (!is.null(benchmark)) {
    # Benchmark the PER-CELL-order point (hat_pc) so the benchmarked point sits on
    # the same per-cell back-transform order as sim_truth. The non-benchmarked
    # point (domains_pred$hat, aggregate order) is left untouched.
    domains_pred$hat_bench <- add_benchmark(x=domains_pred$hat_pc,benchmark_level=benchmark_level,
                                            fwk=fwk,fixed=fixed,benchmark=benchmark,benchmark_type=benchmark_type)
    domains_pred$hat_bench_t <- transform_outcome(domains_pred$hat_bench)$y
    # construct residuals from transformed benchmark estimate
    domains_direct_t <- transform_outcome(collapse:::fmean(sub_domains_direct$outcome,g=sub_domains_direct$domains,w=sub_domains_direct$smp_weights))$y
    domains_direct_t <- data.frame(domains_direct_t = domains_direct_t,domains=names(domains_direct_t))
    domains_direct_t <- dplyr:::inner_join(domains_direct_t,domains_pred,by="domains")
    domains_direct_t_union <- domains_direct_t[domains_direct_t$domains %in% xgb_model$domains,] # limit to domains contains in sample and population
    resid_domains_bench <-  domains_direct_t_union$domains_direct_t - domains_direct_t_union$hat_bench_t
  }

  # Bootstrap
  #_____________________________________________________________________________
  # initialize matrices to store results

  B_results <- matrix(data = NA, nrow = B, ncol = nrow(domains_pred))
  B_results_bench <- matrix(data = NA, nrow = B, ncol = nrow(domains_pred))
  colnames(B_results) <- domains_pred$domains

  # Sample with probability proportional to sample weight if weightsBS = TRUE
  if (weightedBS==T & !is.null(smp_weights)) {
    #pop_subarea_d <- sub_domains_direct$smp_weights/sum(sub_domains_direct$smp_weights)
    pop_subarea_d <- xgb_model$wt_sub_domains/sum(xgb_model$wt_sub_domains)
    #domain_wts <- aggregate(sub_domains_direct$smp_weights,by=list(sub_domains_direct$domains),FUN=sum)$x
    #domain_wts <- collapse:::fsum(sub_domains_direct$smp_weights,g=sub_domains_direct$domains,use.g.names=T)
    domain_wts <- xgb_model$wt_domains
    pop_area_d <- domain_wts/sum(domain_wts)
  } else {
    # Otherwise sample each subarea and area with fixed proobability
    pop_subarea_d <- rep(1/length(resid_sub_domains),length(resid_sub_domains))
    pop_area_d <- rep(1/length(resid_domains),length(resid_domains))
  }

  B_sub <- xgb_model$sub_predictions

  # Standard error of the benchmark-level (e.g. province) direct target, used as
  # the SD of the benchmark-target perturbation, then pre-generate all B
  # perturbations as a named matrix so parallel foreach workers receive a simple
  # object (no rnorm inside workers).
  #
  # Preferred: the caller supplies `benchmark_target_se`, a named vector giving
  # the design-based SE of each benchmark-level direct estimate (keyed by the
  # benchmark_level group value). This should be the SAME estimator that produces
  # the published direct-estimate confidence intervals -- typically the
  # cluster-robust (Taylor) SE -- so the perturbation injects exactly the
  # uncertainty the target is reported with, by construction.
  #
  # Fallback (benchmark_target_se = NULL): the internal Horvitz-Thompson centered
  # variance of the weighted mean (ht_var_weighted_mean), computed from the
  # supplied sample. This is design-based but is computed on the (enumeration-area
  # aggregated) sample handed to xgb(), so it can differ slightly from a
  # cluster-robust SE built on the full microdata; prefer benchmark_target_se when
  # the matching direct SEs are available.
  if (perturb_benchmark && !is.null(benchmark_level)) {
    groups_bm <- unique(as.character(fwk$smp_data[[benchmark_level]]))
    if (!is.null(benchmark_target_se)) {
      nm <- names(benchmark_target_se)
      if (is.null(nm)) stop("perturb_benchmark: benchmark_target_se must be a named vector keyed by benchmark_level group.")
      bm_ht_se <- setNames(rep(NA_real_, length(groups_bm)), groups_bm)
      hit <- intersect(groups_bm, nm)
      bm_ht_se[hit] <- as.numeric(benchmark_target_se[hit])
      if (any(is.na(bm_ht_se))) stop(sprintf(
        "perturb_benchmark: benchmark_target_se is missing %d of %d benchmark groups: %s",
        sum(is.na(bm_ht_se)), length(bm_ht_se),
        paste(names(bm_ht_se)[is.na(bm_ht_se)], collapse = ", ")))
      if (verbose) message(sprintf(
        "perturb_benchmark: supplied direct SE at benchmark level - median %.4f, range [%.4f, %.4f]",
        median(bm_ht_se), min(bm_ht_se), max(bm_ht_se)))
    } else {
      bm_ht_var <- ht_var_weighted_mean(
        y = fwk$smp_data[[fwk$outcome]],
        w = fwk$smp_data[[benchmark_weights]],
        g = fwk$smp_data[[benchmark_level]]
      )
      bm_ht_se <- sqrt(bm_ht_var)
      if (any(is.na(bm_ht_se))) {
        warning("perturb_benchmark: some benchmark groups have n < 2; HT SE set to 0 for those groups.")
        bm_ht_se[is.na(bm_ht_se)] <- 0
      }
      if (verbose) message(sprintf(
        "perturb_benchmark: internal HT SE at benchmark level - median %.4f, range [%.4f, %.4f]",
        median(bm_ht_se), min(bm_ht_se), max(bm_ht_se)
      ))
    }
    # B x n_groups matrix: row j gives the noise to add in bootstrap iteration j.
    # Pre-generating avoids rnorm() inside parallel workers (unpredictable RNG state)
    # and makes the draws reproducible given the outer seed.
    bm_perturbations <- matrix(
      rnorm(B * length(bm_ht_se), mean = 0, sd = rep(bm_ht_se, each = B)),
      nrow = B, ncol = length(bm_ht_se),
      dimnames = list(NULL, names(bm_ht_se))
    )
  } else {
    bm_perturbations <- NULL
  }

  if (bootstrap==T) {

    if (center_residuals==T) {
      resid_sub_domains <- resid_sub_domains-weighted.mean(resid_sub_domains,w=pop_subarea_d)
      resid_domains <- resid_domains-weighted.mean(resid_domains,w=pop_area_d)
      #resid_sub_domains_rescaled <- resid_sub_domains_rescaled-weighted.mean(resid_sub_domains_rescaled,w=pop_subarea_d)
      #resid_domains_rescaled <- resid_domains_rescaled-weighted.mean(resid_domains_rescaled,w=pop_area_d)
      if (!is.null(benchmark)) {
        resid_domains_bench <- resid_domains_bench - weighted.mean(resid_domains_bench,w=pop_area_d)
      }
    }

    if (!is.null(fwk$variance_y)) {
      het_correction <- fwk$smp_data[,fwk$variance_y]^-0.5
    }
    else {
      het_correction <- rep(1,length(fwk$smp_weights_vec))
    }


    clusters <- unique(smp_data[,fwk$domains])
    if (cpus>1) {
      cl <- parallel::makeCluster(cpus)
      doSNOW::registerDoSNOW(cl)
    }
    else {
      foreach::registerDoSEQ()
    }

    # Reproducible parallel RNG. Without this, the foreach %dopar% below draws
    # unseeded random numbers in each worker, so Var_bench and the benchmarked
    # CIs drift run-to-run at the same seed. registerDoRNG assigns each iteration
    # j its own L'Ecuyer-CMRG stream derived from `seed`, independent of how
    # iterations are scheduled across workers, making the bootstrap byte-for-byte
    # reproducible. Applies to both the doSNOW (parallel) and doSEQ (sequential)
    # backends and does not affect the point-estimate RNG (set.seed elsewhere).
    doRNG::registerDoRNG(seed)

    cat("Beginning bootstrap \n")
    #for (j in 1:B){

    pb <- txtProgressBar(max = B, style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)

    #browser()
    B_results_list <- foreach(j = 1:B,
                              .combine = rbind,
                              .errorhandling = "stop",
                              .packages = c("xgboost"),
                              .options.snow = list(progress = progress)) %dopar% {


                                #displayevery = max(round(B/10,0),1)

                                #if (j %% displayevery==0) {
                                #cat(paste0("replication ",j," of ",B,"\n"))
                                #}

                                if (bootstrap_type=="residual") {
                                  # First implement a standard residual bootstrap
                                  # randomly sample residuals USING THE WEIGHTS CALCULATED ABOVE if weighted_BS==TRUE
                                  sub_domains_draw <- resid_sub_domains[sample(1:length(resid_sub_domains),nrow(B_sub), prob=pop_subarea_d,
                                                                               replace = TRUE)]
                                  # This undoes rescaling of subarea residuals if variance_y is specified
                                  #if (!is.null(variance_y)) {
                                  #sub_domains_draw <- sub_domains_draw*(B_sub[,fwk$pop_weights]^-0.5)
                                  #}

                                  B_sub$sim_t <- as.numeric(B_sub$hat_t) + sub_domains_draw
                                  # aggregate to area
                                  B_domains <- collapse:::fmean(x=B_sub$sim_t,g=B_sub[,fwk$domains],w=B_sub[,fwk$pop_weights])

                                  B_domains <- data.frame("domains" = names(B_domains),"sim_t" = B_domains)
                                  if (!is.null(benchmark_level)) {
                                    B_domains[,benchmark_level]<-collapse:::ffirst(x=B_sub[,benchmark_level],g=B_sub[,fwk$domains])
                                  }
                                  B_domains$hat_t <- collapse:::fmean(x=B_sub$hat_t,g=B_sub$domains,w=B_sub$wts)

                                  #B_domains <- data.frame("domains" = unique(B_sub$domains),"sim" = B_domains)
                                  area_draw <- data.frame(domains=B_domains$domains,area_draw=resid_domains[sample(1:length(resid_domains),
                                                                                                                   nrow(B_domains), prob=pop_area_d,
                                                                                                                   replace = TRUE)])
                                  B_domains <- left_join(B_domains,area_draw,by="domains")

                                  # Add bootstrapped area residual and back transform
                                  B_domains$sim_t_plus_area <- B_domains$sim_t +B_domains$area_draw
                                  B_domains$sim_plus_area <- back_transform_outcome(B_domains$sim_t_plus_area)

                                  # [ORDER] Per-cell replicate, so the interval quantiles follow the
                                  # same back-transform-then-aggregate order as the point estimate
                                  # (hat_pc). Mirrors point_estim_xgb lines ~1144-1149 and the
                                  # benchmark branch below. The area draw is REUSED, not redrawn, so
                                  # this consumes no additional RNG and the bootstrap stream is
                                  # unchanged; with transformation = "no" back_transform_outcome is
                                  # the identity and sim_plus_area_pc == sim_plus_area exactly.
                                  sub_area_draw_pc <- area_draw$area_draw[match(B_sub[,fwk$domains],
                                                                                area_draw$domains)]
                                  B_domains$sim_plus_area_pc <- collapse:::fmean(
                                    x = back_transform_outcome(B_sub$sim_t + sub_area_draw_pc),
                                    g = B_sub[,fwk$domains], w = B_sub[,fwk$pop_weights])



                                  if (!is.null(benchmark)) {
                                    # Benchmarking requires a different bootstrapping procedure because the residuals are no longer independent.
                                    # We implement the same type of approach used in parametric bootstraps
                                    # 1. Redraw areas from benchmarked area distribution and calculate simulated "truth" mean for each domain

                                    # Draw from resid_domains (model-based area residuals) rather than resid_domains_bench.
                                    # resid_domains_bench is (direct estimate - benchmarked model prediction) and is dominated
                                    # by sampling error in the direct estimate, not by the area random effect. The bootstrap
                                    # population should reflect the model's data-generating process, which uses resid_domains.
                                    area_draws_bench <- data.frame(domains=B_domains$domains,area_draw=resid_domains[sample(1:length(resid_domains),
                                                                                                                            nrow(B_domains), prob=pop_area_d,
                                                                                                                            replace = TRUE)])

                                    colnames(area_draws_bench)[1]=fwk$domains
                                    B_sub <- dplyr:::left_join(B_sub,area_draws_bench,by=fwk$domains)
                                    B_sub$sim_t_plus_area <- B_sub$sim_t+B_sub$area_draw
                                    B_sub$sim_plus_area <- back_transform_outcome(B_sub$sim_t_plus_area)
                                    B_sub$area_draw <- NULL

                                    B_domains$sim_truth <- collapse:::fmean(B_sub$sim_plus_area,g=B_sub$domains,w=B_sub$wts)

                                    # 2. Extract sample observations from new population
                                    covariates_no_popwt <- setdiff(fwk$covariates, fwk$pop_weights)
                                    B_sample <- dplyr:::inner_join(smp_data[, c(fwk$sub_domains, covariates_no_popwt, fwk$smp_weights, fwk$variance_y)]
                                                                   , B_sub, by=fwk$sub_domains)
                                    #B_sample$sim_plus_area <- back_transform_outcome(B_sample$sim_t_plus_area)

                                    #2.5 Truncate new dependent variable if necessary
                                    if ((!is.null(list(...)$objective) && list(...)$objective=="reg:logistic") & transformation=="no") {
                                      B_sample$sim_t_plus_area <- pmax(0, pmin(1, B_sample$sim_t_plus_area))
                                    }
                                    else if (!is.null(list(...)$objective) && list(...)$objective=="reg:gamma" && bootstrap_type=="residual") {
                                      B_sample$sim_t_plus_area <- pmax(1e-5, B_sample$sim_t_plus_area)
                                    }
                                    else if (transformation=="arcsin") {
                                      # Clamp to [0, 1] to prevent NaN from arcsin transform of out-of-range values
                                      B_sample$sim_plus_area <- pmax(0, pmin(1, B_sample$sim_plus_area))
                                    }
                                    else if (transformation=="log") {
                                      # Ensure positive values for log transform
                                      B_sample$sim_plus_area <- pmax(1e-10, B_sample$sim_plus_area)
                                    }
                                    else if (transformation=="log.shift") {
                                      # Ensure y + ls_lambda > 0 for log(y + ls_lambda)
                                      B_sample$sim_plus_area <- pmax(-ls_lambda + 1e-10, B_sample$sim_plus_area)
                                    }

                                    # Guard: error if any NaN/Inf in the bootstrap dependent variable before re-estimation
                                    boot_y <- if (transformation == "arcsin") B_sample$sim_plus_area else B_sample$sim_t_plus_area
                                    bad_vals <- is.nan(boot_y) | is.infinite(boot_y)
                                    if (any(bad_vals)) {
                                      stop(sprintf("Bootstrap iteration produced %d NaN/Inf values in simulated outcome. This may indicate that benchmarking is producing out-of-range estimates.", sum(bad_vals)))
                                    }

                                    X_smp_boot <- B_sample[,fwk$covariates]
                                    if (!is.null(fwk$variance_y)) {
                                      boot_het <- B_sample[,fwk$variance_y]^-0.5
                                    } else {
                                      boot_het <- 1
                                    }
                                    boot_weights <- B_sample[,fwk$smp_weights] * boot_het
                                    # Rescale weights within each domain so they sum to the domain sample size
                                    if (rescale_weights) {
                                      boot_domain_sum_wts <- ave(boot_weights, B_sample[,fwk$domains], FUN = sum)
                                      boot_domain_n <- ave(boot_weights, B_sample[,fwk$domains], FUN = length)
                                      boot_weights <- boot_weights * boot_domain_n / boot_domain_sum_wts
                                    }
                                    boot_weights <- boot_weights / mean(boot_weights)

                                    # 3. Generate predictions using new sample values in arcsin space
                                    predictions  <- point_estim_xgb(params=params,smp_X=X_smp_boot,
                                                                    smp_Y=B_sample$sim_plus_area,
                                                                    smp_weight= boot_weights,
                                                                    pop_X=X_pop_xgb,
                                                                    sub_domains=sub_domains,
                                                                    nrounds=nrounds,
                                                                    fwk=fwk,
                                                                    transform_outcome=transform_outcome,
                                                                    back_transform_outcome=back_transform_outcome,
                                                                    L=L,variance_y=NULL)$predictions



                                    # benchmark predictions to the collapsed simulated sample in probability space
                                    # First calculate weighted mean of sample to benchmark, for each benchmark_level
                                    B_sample_bm <- data.frame(sample_bm = collapse:::fmean(x=B_sample$sim_plus_area,g=B_sample[,benchmark_level],w=B_sample[,benchmark_weights]))
                                    B_sample_bm[,benchmark_level] <- rownames(B_sample_bm)

                                    # (benchmark-target perturbation is applied post-loop in the main process)

                                    #B_domains <- dplyr:::left_join(B_domains,domains_pred[,c("domains","hat_bench","weight")],by="domains")
                                    #B_domains$mean_sim_plus_area <- collapse:::fmean(B_domains$sim_plus_area,g=B_domains[,benchmark_level],w=B_domains$weight,TRA=1)
                                    #B_domains <- dplyr:::left_join(B_domains,B_sample_bm,by=benchmark_level)

                                    bm_vec <- setNames(B_sample_bm$sample_bm, B_sample_bm[,benchmark_level])
                                    # Benchmark the PER-CELL-order replicate point (hat_pc) so each
                                    # benchmarked replicate matches sim_truth's order; Var_bench then
                                    # differences like-order quantities. sim_truth (per-cell) is unchanged.
                                    B_domains$sim_bench <- add_benchmark(predictions$hat_pc,benchmark_level=benchmark_level,benchmark=bm_vec,fwk=fwk,
                                                                         fixed=fixed,benchmark_type=benchmark_type)
                                  } # close additional benchmarking bootstrap code
                                } # close residual bootstrap code to produce B_domains_sim

                                # case resampling bootstrap
                                else if (bootstrap_type=="case") {
                                  # Resample of clusters with replacement
                                  #browser()
                                  boot_clusters <- sample(clusters, replace = TRUE)
                                  # Get all observations from sampled clusters
                                  boot_data <- do.call(rbind, lapply(boot_clusters, function(clust) {
                                    smp_data[smp_data[,fwk$domains] == clust, ]}))
                                  X_smp_boot <- boot_data[fwk$covariates]
                                  # Rescale weights within each domain so they sum to the domain sample size
                                  case_wts <- boot_data[,smp_weights]
                                  if (rescale_weights) {
                                    case_domain_sum_wts <- ave(case_wts, boot_data[,fwk$domains], FUN = sum)
                                    case_domain_n <- ave(case_wts, boot_data[,fwk$domains], FUN = length)
                                    case_wts <- case_wts * case_domain_n / case_domain_sum_wts
                                  }
                                  #Estimate model
                                  dtrain <- xgboost::xgb.DMatrix(
                                    data = as.matrix(X_smp_boot),
                                    label = transform_outcome(boot_data[,fwk$outcome])$y,
                                    weight = (case_wts/mean(case_wts))
                                  )
                                  xgb_fit <- xgboost::xgb.train(
                                    data               = dtrain,
                                    params             = params,
                                    nrounds            = nrounds,
                                    verbose            = 0
                                  )
                                  # Obtain predictions
                                  B_sub$sim <-
                                    back_transform_outcome(predict(xgb_fit, as.matrix(X_pop_xgb))
                                    )
                                  B_domains <- collapse:::fmean(x=B_sub$sim,g=B_sub$domains,w=B_sub$wts)
                                  B_domains <- data.frame("domains" = names(B_domains),"sim" = B_domains)
                                  B_domains$sim_plus_area <- B_domains$sim
                                  # [ORDER] this branch already back-transforms per cell before
                                  # aggregating, so the two orderings coincide here
                                  B_domains$sim_plus_area_pc <- B_domains$sim
                                  B_domains$sim_bench <- NULL
                                  B_domains$sim_truth <- NULL
                                }


                                # Return a single-row data frame with list-columns
                                data.frame(
                                  B_results = I(list(B_domains$sim_plus_area)),      # I() prevents unlisting
                                  B_results_pc = I(list(B_domains$sim_plus_area_pc)),  # [ORDER]
                                  B_results_bench = I(list(B_domains$sim_bench)),
                                  B_results_truth = I(list(B_domains$sim_truth)),
                                  B_domains = I(list(B_domains$domains)),
                                  stringsAsFactors = FALSE)

                              } # close bootstrap loop

    close(pb)

    if (cpus>1) {
      parallel::stopCluster(cl)
    }

    # Extract and combine from list

    B_results <- do.call(rbind, B_results_list$B_results)
    B_results_pc <- do.call(rbind, B_results_list$B_results_pc)   # [ORDER]
    B_results_bench <- do.call(rbind, B_results_list$B_results_bench)
    B_results_truth <- do.call(rbind, B_results_list$B_results_truth)
    B_results_domains <-  B_results_list$B_domains[[1]]


    #browser()
    # This procedure resorts the data when using multiple cores, so we will undo the resorting
    original_domain_order <- unique(fwk$pop_domains_vec)
    row_reorder <- match(original_domain_order, B_results_domains)
    B_results <- B_results[,row_reorder]
    B_results_pc <- B_results_pc[,row_reorder]   # [ORDER]
    if (!is.null(benchmark)) {
      B_results_bench <- B_results_bench[,row_reorder]
      B_results_truth <- B_results_truth[,row_reorder]
    }
    B_results_domains <- B_results_domains[row_reorder]
    #tail(domains_pred)
    #tail(B_results_domains)

    # Post-loop benchmark-target perturbation (runs in main process — no parallel
    # scoping issues). For each state s and each bootstrap iteration j, adds the
    # pre-generated noise δ[j,s] to the benchmarked replicate of every ward in s.
    # This propagates survey-sampling uncertainty of the state direct estimate into
    # the benchmarked CIs. Equivalent to within-loop perturbation of B_sample_bm
    # under the linear approximation that the benchmark adjusts wards additively.
    if (!is.null(bm_perturbations) && !is.null(benchmark)) {
      # Save unperturbed per-ward Var_bench BEFORE adding noise / clipping, so it
      # is not censored by the [0,1] clip applied to perturbed replicates below.
      # Useful for ward-vs-ward comparisons where the state-level shift cancels.
      var_bench_unperturbed <- apply(B_results_bench - B_results_truth, 2, var)

      # Build ward → benchmark-group (state) lookup from B_sub (has both columns)
      grp_raw   <- collapse:::ffirst(x = B_sub[[benchmark_level]], g = B_sub[[fwk$domains]])
      dom_to_grp <- setNames(as.character(grp_raw), names(grp_raw))
      grp_per_col <- dom_to_grp[as.character(B_results_domains)]  # one entry per column

      if (verbose) {
        var_before <- mean(apply(B_results_bench, 2, var), na.rm = TRUE)
        n_matched  <- sum(!is.na(grp_per_col) & grp_per_col %in% colnames(bm_perturbations))
        message(sprintf(
          "perturb_benchmark DIAG: %d/%d columns matched to a benchmark group; mean col-var BEFORE = %.6f",
          n_matched, ncol(B_results_bench), var_before
        ))
      }

      for (grp in unique(grp_per_col[!is.na(grp_per_col)])) {
        if (!grp %in% colnames(bm_perturbations)) next
        cols <- which(grp_per_col == grp)
        # Perturb the benchmark target on the model's transformed (unbounded) scale
        # and back-transform, so perturbed draws stay in range WITHOUT a clamp. The
        # rate-scale target SE (bm_perturbations) is mapped to a transformed-scale
        # increment via the transform's local derivative at the group mean (delta
        # method), preserving the intended benchmark-uncertainty magnitude. The
        # previous behaviour added the perturbation on the rate scale and hard-clamped
        # arcsin draws to [0, 1]; near the 0/1 boundary that clamp could CONTRACT
        # interval widths (occasionally below the unperturbed width), an artefact that
        # perturbing on the transformed scale removes. bm_perturbations[, grp] is
        # length-B and is broadcast across all cols.
        if (transformation == "no") {
          B_results_bench[, cols] <- B_results_bench[, cols] + bm_perturbations[, grp]
        } else {
          lo   <- 1e-9
          hi   <- if (transformation == "arcsin") 1 - 1e-9 else Inf
          pbar <- min(max(mean(B_results_bench[, cols], na.rm = TRUE), lo), hi)
          eps  <- 1e-4
          p_hi <- min(pbar + eps, hi); p_lo <- max(pbar - eps, lo)
          dzdp <- (transform_outcome(p_hi)$y - transform_outcome(p_lo)$y) / (p_hi - p_lo)
          tb   <- transform_outcome(pmin(pmax(B_results_bench[, cols], lo), hi))$y
          tb   <- tb + bm_perturbations[, grp] * as.numeric(dzdp)
          B_results_bench[, cols] <- back_transform_outcome(tb)
        }
      }

      if (verbose) {
        var_after <- mean(apply(B_results_bench, 2, var), na.rm = TRUE)
        message(sprintf(
          "perturb_benchmark DIAG: mean col-var AFTER  = %.6f  (delta = %+.6f, expected approx %.6f)",
          var_after, var_after - var_before, mean(bm_ht_se^2)
        ))
      }
    } else {
      var_bench_unperturbed <- NULL
      if (verbose) message(sprintf(
        "perturb_benchmark DIAG: SKIPPED — bm_perturbations is %s, benchmark is %s",
        ifelse(is.null(bm_perturbations), "NULL", "set"),
        ifelse(is.null(benchmark), "NULL", "set")
      ))
    }

    # Prepare results
    #_____________________________________________________________________________
    #results <- NULL
    results <- data.frame(domains = B_results_domains)
    results$Mean_boot <- NA
    results$Mean_boot_agg <- NA   # [ORDER]
    results$Lower_boot <- NA
    results$Upper_boot <- NA
    results$Var_boot <- NA
    results$Mean_boot_bench <- NA
    results$Lower_bench <- NA
    results$Upper_bench <- NA

    if (!is.null(benchmark)) {
      results <- left_join (results, domains_pred[, c("hat","hat_pc","hat_bench","domains")],by="domains")
    } else {
      results <- left_join (results, domains_pred[, c("hat","hat_pc","domains")],by="domains")
    }

    colnames(results)[1] <- "Domain"
    #browser()

    for (l in 1:ncol(B_results)){

      # [ORDER] quantiles are taken from the per-cell replicate so the interval and the
      # reported point estimate (hat_pc) are on the same back-transform-then-aggregate
      # order. B_results (aggregate-then-back-transform) is retained for Mean_boot_agg.
      temp <- B_results_pc[,l]
      temp_agg <- B_results[,l]

      # if (transformation=="arcsin"){
      #   temp <- ifelse(temp>asin(1), asin(1), temp)
      #   temp <- ifelse(temp<asin(0), asin(0), temp)
      #   temp <- sin(temp)^2
      # }
      # if (transformation=="log"){
      #   temp <- exp(temp)
      # }
      results$Mean_boot[l] <- mean(temp)
      results$Mean_boot_agg[l] <- mean(temp_agg)   # [ORDER] old ordering, for reconciliation
      results$Lower_boot[l] <- quantile(temp, probs = (1-conf_level)/2)
      results$Upper_boot[l] <- quantile(temp, probs = 1-(1-conf_level)/2)
      results$Var_boot[l] <- var(temp)
      if (!is.null(benchmark)) {
        resid <- B_results_bench[,l]-B_results_truth[,l]
        resid <- resid - mean(resid)
        # Report var(resid) after centering rather than mean(resid^2). The squared bias term
        # in mean(resid^2) is dominated by a structural offset between the simulated truth and
        # the benchmarked estimator (the truth aggregates raw model predictions plus residuals,
        # while the estimator applies a multiplicative benchmark adjustment), not by real
        # estimator bias. var(resid) gives a variance summary consistent with the CI quantiles below.
        results$Var_bench[l] <- var(resid)
        temp_bench <- B_results_bench[,l]
        results$Mean_boot_bench[l] <- mean(temp_bench)
        results$Lower_bench[l] <- results$hat_bench[l] + quantile(resid, probs = (1-conf_level)/2)
        results$Upper_bench[l] <- results$hat_bench[l] + quantile(resid, probs = 1-(1-conf_level)/2)
      } # Close benchmark
    } # close loop over areas

  } # close if bootstrap branch

  #browser()
  sub_domains_direct$outcome_t <- transform_outcome(sub_domains_direct$outcome)$y
  sub_domains_direct$hat <- NA
  in_pop <- sub_domains_direct[,sub_domains] %in% fwk$X_pop[,sub_domains]
  sub_domains_direct$hat[in_pop] <- back_transform_outcome(sub_domains_direct$outcome_t[in_pop] - xgb_model$resid_total)

  #sub_domains_direct$hat <- back_transform_outcome(sub_domains_direct$outcome_t-xgb_model$resid_total)

  if (bootstrap==T) {
    # [ORDER] Mean is hat_pc (back-transform per cell, then aggregate), which is the
    # coherent order and matches what the benchmarked path has always used. hat
    # (aggregate on the transformed scale, then back-transform) is retained as
    # Mean_agg so results from earlier versions can be reconciled. With
    # transformation = "no" the two are identical.
    # The CI recentring uses hat_pc and Mean_boot, both per-cell, so point estimate
    # and interval are on the same order.
    result <- list(
      ind = data.frame(Domain = results$Domain, Mean = results$hat_pc,
                       Mean_agg = results$hat),
      var = data.frame(Domain = results$Domain, Mean = results$Var_boot),
      CI  = data.frame(Domain = results$Domain,
                       Lower = results$Lower_boot+results$hat_pc-results$Mean_boot,
                       Upper = results$Upper_boot+results$hat_pc-results$Mean_boot),
      yhat=data.frame(sub_domains_direct[,c(sub_domains,"hat")]),
      model = xgb_fit,
      smp_data =  smp_data,
      out_call = out_call,
      transformation = transformation,
      framework = fwk
    )


    if (boot_estimates==T)  {
      result$ind[,2] <- results["Mean_boot"]
      result$CI$Lower <- result$CI$Lower+results["Mean_boot"]-results["hat"]
      result$CI$Upper <- result$CI$Upper+results["Mean_boot"]-results["hat"]
    } # close use boot estimates
    colnames(result$var)[2]="Mean"
    colnames(result$ind)[2]="Mean"
  } # close bootstrap==T
  else {
    # Bootstrap not selected
    result <- list(
      # [ORDER] see the bootstrap branch above: Mean is the per-cell order, Mean_agg
      # preserves the previous aggregate-then-back-transform value.
      ind = data.frame(Domain = domains_pred[,"domains"], Mean = domains_pred[,"hat_pc"],
                       Mean_agg = domains_pred[,"hat"]),
      yhat=data.frame(sub_domains_direct[,c(sub_domains,"hat")]),
      model = xgb_fit,
      smp_data =  smp_data,
      out_call = out_call,
      transformation = transformation,
      framework = fwk
    )
  } # close bootstrap not selected

  if (!is.null(benchmark)) {
    if (bootstrap==T) {
      result$ind$Mean_bench <- results$hat_bench
      result$CI$Lower_bench <- results$Lower_bench
      result$CI$Upper_bench <- results$Upper_bench
      if (boot_estimates==T) {
        result$ind$Mean_bench <- results$Mean_boot_bench
        result$CI$Lower_bench <- results$Lower_bench-results["hat_bench"]-results["Mean_boot_bench"]
        result$CI$Upper_bench <- results$Upper_bench-results["hat_bench"]-results["Mean_boot_bench"]
      }
      result$var$Var_bench <- results$Var_bench
      if (!is.null(result$var)) {
        colnames(result$var)[ncol(result$var)] <- "Var_bench"
      }
    } else {
      # Add benchmarked estimates even without bootstrap
      result$ind$Mean_bench <- domains_pred$hat_bench
    }
    # Per-ward Var_bench computed BEFORE benchmark-target perturbation/clipping.
    # Use this for ward-vs-ward comparisons (the additive state-level perturbation
    # cancels in within-state contrasts). Stored in its own slot rather than in
    # result$var so it does not break downstream CV/MSE machinery that expects
    # result$var and result$ind to have matching column counts.
    # NULL when perturb_benchmark=FALSE.
    if (bootstrap == TRUE && !is.null(var_bench_unperturbed)) {
      result$Var_bench_unperturbed <- data.frame(
        Domain = result$var$Domain,
        Var_bench_unperturbed = var_bench_unperturbed[
          match(result$var$Domain, B_results_domains)
        ]
      )
    } else {
      result$Var_bench_unperturbed <- NULL
    }
  } # close benchmark

  # Truncation — applied whenever bootstrap CIs exist
  if (bootstrap==T) {
    if (transformation=="arcsin" | (!is.null(list(...)$objective) && list(...)$objective=="reg:logistic")) {
      result$CI$Lower <- pmax(result$CI$Lower, 0)
      result$CI$Upper <- pmin(result$CI$Upper, 1)
      if (!is.null(benchmark)) {
        result$CI$Lower_bench <- pmax(result$CI$Lower_bench, 0)
        result$CI$Upper_bench <- pmin(result$CI$Upper_bench, 1)
      }
    }
    else if (transformation=="poisson" || (!is.null(list(...)$objective) && list(...)$objective=="reg:gamma" && bootstrap_type=="residual")) {
      result$CI$Lower <- pmax(result$CI$Lower, 0)
      if (!is.null(benchmark)) {
        result$CI$Lower_bench <- pmax(result$CI$Lower_bench, 0)
      }
    }
  }


  class(result) <- c("xgb","povmap")
  return(result)

} # close xgb function


point_estim_xgb <- function(params, smp_X, smp_Y, smp_weight, pop_X, sub_domains,nrounds,fwk,transform_outcome,back_transform_outcome,L,variance_y) {

  smp_Y_t <- transform_outcome(smp_Y)$y

  if (!is.null(variance_y)) {
    het_correction <- fwk$smp_data[,variance_y]^-0.5
  }
  else {
    het_correction <- rep(1,length(smp_weight))
  }

  dtrain <- xgboost:::xgb.DMatrix(
    data = as.matrix(smp_X),
    label = smp_Y_t,
    weight = smp_weight*het_correction/mean(smp_weight*het_correction)
  )

  xgb_fit <- xgboost:::xgb.train(
    data               = dtrain,
    params             = params,
    nrounds            = nrounds,
    verbose            = 0

  )



  #Generate predictions at sub-area level in transformed space
  sub_pred_t <- data.frame(
    fwk$X_pop[,sub_domains],
    fwk$X_pop[,fwk$pop_weights],
    fwk$X_pop[,fwk$domains],
    predict(xgb_fit, as.matrix(pop_X))
  )
  if (!is.null(fwk$benchmark_level)) {
    sub_pred_t <- data.frame(
      fwk$X_pop[,sub_domains],
      fwk$X_pop[,fwk$pop_weights],
      fwk$X_pop[,fwk$domains],
      fwk$X_pop[,fwk$benchmark_level],
      predict(xgb_fit, as.matrix(pop_X))
    )
    colnames(sub_pred_t) <- c(fwk$sub_domains, fwk$pop_weights, fwk$domains, fwk$benchmark_level, "hat_t")
  } else {
    colnames(sub_pred_t) <- c(fwk$sub_domains, fwk$pop_weights, fwk$domains, "hat_t")
  }

  # calculate sub-area residuals
  sample <- inner_join(fwk$smp_data[,c(fwk$sub_domains,fwk$outcome,fwk$smp_weights,variance_y)],sub_pred_t,by=fwk$sub_domains)
  outcome_t <- paste0(fwk$outcome,"_t")
  sample[,outcome_t] <- transform_outcome(sample[,fwk$outcome])$y
  resid_total <- sample[,outcome_t]-sample$hat_t

  # rescale residuals if variance_y is provided
  #  if (!is.null(fwk$variance_y)) {
  #resid_total_rescaled <- resid_total/(sample[,fwk$variance_y]^0.5)
  #} else {
  #  resid_total_rescaled <- resid_total
  #  }

  wt_sub_domains <- sample[,fwk$smp_weights]
  sample_domains <- data.frame(direct=collapse:::fmean(sample[,outcome_t],g=sample[,fwk$domains],w=wt_sub_domains),
                               hat=collapse:::fmean(sample$hat_t,g=sample[,fwk$domains],w=wt_sub_domains),
                               wts= collapse:::fsum(wt_sub_domains,g=sample[,fwk$domains]),
                               resid_domains=collapse:::fmean(resid_total,g=sample[,fwk$domains],w=wt_sub_domains))

  sample_domains <- data.frame(sample_domains,rownames(sample_domains))
  colnames(sample_domains)[ncol(sample_domains)] <- fwk$domains
  sample <- left_join(sample,sample_domains,by=fwk$domains)

  resid_domains <- sample_domains$resid_domains
  resid_sub_domains <- resid_total -sample$resid_domains

  #resid_domains_rescaled <- sample_domains$resid_domains_rescaled
  #resid_sub_domains_rescaled <- resid_total_rescaled - sample$resid_domains_rescaled

  mean_sim    <- 0   # aggregate-then-back-transform (non-benchmarked point order)
  mean_sim_pc <- 0   # per-cell-then-aggregate (benchmarked point order; matches sim_truth)



  for (l in 1:L) {
    sub_pred_t$sim_t <- as.numeric(sub_pred_t$hat_t) + resid_sub_domains[sample(1:length(resid_sub_domains),
                                                                                nrow(sub_pred_t), prob=sample[,fwk$smp_weights],
                                                                                replace = TRUE)]
    B_domains <- collapse:::fmean(x=sub_pred_t$sim_t,g=sub_pred_t[,fwk$domains],w=sub_pred_t[,fwk$pop_weights])
    B_domains <- data.frame("domains" = names(B_domains),"sim_t" = B_domains)
    B_domains <-data.frame(B_domains,wts=collapse:::fsum(x=sub_pred_t[,fwk$pop_weights],g=sub_pred_t[,fwk$domains]))
    area_draws <- data.frame(domains=B_domains$domains,area_draw=resid_domains[sample(1:length(resid_domains),
                                                                                      nrow(B_domains), prob=sample_domains$wts,
                                                                                      replace = TRUE)])

    B_domains <- left_join(B_domains,area_draws,by="domains")

    # Aggregate-then-back-transform: the NON-benchmarked point order (unchanged).
    B_domains$sim_t_plus_area <- B_domains$sim_t +B_domains$area_draw
    B_domains$sim_plus_area <- back_transform_outcome(B_domains$sim_t_plus_area)
    mean_sim <- (mean_sim*(l-1)+B_domains$sim_plus_area)/l

    # Per-cell-then-aggregate: the BENCHMARKED point order, constructed exactly as
    # sim_truth in the bootstrap -- broadcast each area's draw to its cells, add on
    # the transformed scale, back-transform EACH CELL, then population-weight
    # aggregate. No extra RNG is drawn here, so the aggregate path (mean_sim) and
    # the whole RNG stream stay byte-for-byte unchanged.
    sub_area_draw <- area_draws$area_draw[match(sub_pred_t[,fwk$domains], area_draws$domains)]
    sim_pa_cell   <- back_transform_outcome(sub_pred_t$sim_t + sub_area_draw)
    pc <- collapse:::fmean(x=sim_pa_cell,g=sub_pred_t[,fwk$domains],w=sub_pred_t[,fwk$pop_weights])
    mean_sim_pc <- (mean_sim_pc*(l-1) + as.numeric(pc[as.character(B_domains$domains)]))/l
  } # end back-transformation loop to calculate mean of simulations through back-transformation

  # generate point estimates
  MC_results <-  data.frame(B_domains["domains"],
                            hat=mean_sim,
                            hat_pc=mean_sim_pc)
  sub_pred_t$sim_t <- NULL


  return(list(predictions=MC_results,sub_predictions = sub_pred_t, model=xgb_fit,resid_sub_domains=resid_sub_domains,
              resid_domains=resid_domains, resid_total=resid_total,wt_sub_domains=wt_sub_domains,
              wt_domains=sample_domains$wts,domains=rownames(sample_domains)))
}




add_benchmark <- function(x, benchmark_level,fwk,fixed,benchmark,benchmark_type) {
  point_estim <- NULL
  point_estim$ind <- data.frame(Mean = x)

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
  bench <- point_estim$ind$Mean_bench
  return(bench)
}