#' Cross-validation for extreme gradient boosting SAE models using population data
#'
#' The function \code{xgb_cv} performs k-fold cross-validation for extreme gradient
#' boosting SAE models, holding out entire domains from the sample while retaining
#' the full population data for prediction. This provides a more realistic assessment
#' of out-of-sample performance than sample-only cross-validation, since it mirrors
#' the actual SAE workflow where population covariates are always available.
#'
#' @param fixed a two-sided linear formula object describing the
#' fixed-effects part of the model with the dependent variable on the left
#' of a ~ operator and the explanatory variables on the right, separated
#' by + operators.
#' @param smp_data a data frame containing the sample data.
#' @param smp_weights a character string containing the name of the variable that
#' indicates weights in \code{smp_data}. Defaults to \code{NULL}.
#' @param pop_data a data frame containing the population/census data.
#' @param pop_weights a character string containing the name of the variable that
#' indicates population weights in \code{pop_data}. Defaults to \code{NULL}.
#' @param domains a character string containing the name of a variable
#' that indicates domains in \code{smp_data} and \code{pop_data}.
#' @param sub_domains character string specifying the variable name that denotes
#' sub-domains within the dataset.
#' @param transformation a character string. Transformation types for the dependent
#' variable: "no", "log", "arcsin", "log.shift", or "logit". Defaults to \code{"no"}.
#' @param folds number of cross-validation folds. Defaults to 10.
#' @param L number of Monte-Carlo simulations for back-transformation. Defaults to 100.
#' @param nrounds maximum number of boosting iterations. Defaults to 100.
#' @param max_depth maximum depth of a tree. Defaults to 4.
#' @param colsample_bytree subsample ratio of columns per tree. Defaults to 0.6.
#' @param colsample_bylevel subsample ratio of columns per level. Defaults to 0.6.
#' @param colsample_bynode subsample ratio of columns per node. Defaults to 0.6.
#' @param subsample subsample ratio of training instances. Defaults to 0.6.
#' @param min_child_weight minimum sum of instance weight in a child. Defaults to 1.
#' @param eta step size shrinkage. Defaults to 0.3.
#' @param gamma minimum loss reduction for partitioning. Defaults to 0.
#' @param max_delta_step maximum delta step. Defaults to 0.
#' @param lambda L2 regularization term. Defaults to 1.
#' @param alpha L1 regularization term. Defaults to 0.
#' @param na.rm if TRUE, observations with NA values are removed. Defaults to FALSE.
#' @param seed integer seed for reproducibility. Defaults to 123.
#' @param rescale_weights if TRUE, rescales sample weights within each domain so
#' they sum to the domain sample size. Defaults to TRUE.
#' @param variance_y the name of the variable containing the variance of the outcome.
#' Defaults to NULL.
#' @param cpus number of cores to parallelize across folds. Defaults to 1 (sequential).
#' @param verbose display progress. Defaults to TRUE.
#' @param plot_filename optional file path to save a ggplot2 scatter plot of
#' predicted vs. direct estimates (e.g. "cv_plot.png"). If NULL, no plot is
#' generated. Defaults to NULL.
#' @param ... additional parameters passed to \code{xgb}.
#'
#' @return A list containing (all fit statistics are scored against the direct
#'   estimate, a noisy target, not against a known truth):
#' \describe{
#'   \item{r2_cv}{Out-of-sample R-squared, the proportion of variance explained,
#'     exactly 1 - SSE/SST with matched (sum) denominators. Penalizes bias and
#'     scale compression and can be negative. Denominator is the variation of the
#'     noisy direct estimate.}
#'   \item{cor_cv}{Out-of-sample Pearson correlation between the held-out
#'     prediction and the direct estimate. Invariant to a linear rescaling of the
#'     prediction, so it is not driven down by compression toward the mean.}
#'   \item{cor2_cv}{Square of \code{cor_cv}, the squared Pearson correlation.
#'     Distinct from \code{r2_cv}; the two coincide only for in-sample OLS fits.}
#'   \item{mse_cv}{Out-of-sample MSE on the back-transformed scale}
#'   \item{mae_cv}{Out-of-sample MAE on the back-transformed scale}
#'   \item{rank_cor_cv}{Out-of-sample Spearman rank correlation}
#'   \item{skewness_domains}{Skewness of domain-level outcomes}
#'   \item{domain_results}{Data frame with domain-level predictions and direct estimates}
#'   \item{plot}{ggplot2 object (only if \code{plot_filename} is specified)}
#' }
#'
#' @export

xgb_cv <- function(fixed,
                   smp_data,
                   smp_weights = NULL,
                   pop_data,
                   pop_weights = NULL,
                   domains,
                   sub_domains,
                   transformation = "no",
                   folds = 10,
                   L = 100,
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
                   rescale_weights = TRUE,
                   variance_y = NULL,
                   cpus = 1,
                   verbose = TRUE,
                   plot_filename = NULL,
                   ...) {

  # Pinned-version guard (see BUILD_PIN_xgb.txt) -- fail before any model fit
  .assert_xgb_version()

  set.seed(seed)

  # Extract outcome variable name
  outcome <- all.vars(fixed[[2]])

  # Compute direct estimates from the FULL sample for each domain
  # These serve as the "truth" to compare against
  if (!is.null(smp_weights)) {
    wts_full <- smp_data[, smp_weights]
  } else {
    wts_full <- rep(1, nrow(smp_data))
  }
  domain_vec_full <- smp_data[, domains]
  outcome_full <- smp_data[, outcome]

  direct_df <- data.frame(y = outcome_full, d = domain_vec_full, w = wts_full)
  direct_means <- sapply(split(direct_df, direct_df$d),
                         function(g) weighted.mean(g$y, w = g$w))
  direct_estimates <- data.frame(Domain = names(direct_means),
                                  Direct = as.numeric(direct_means),
                                  stringsAsFactors = FALSE)

  # Assign domains to folds
  unique_domains <- unique(smp_data[, domains])
  fold_assignment <- data.frame(Domain = unique_domains,
                                 fold = sample(rep(1:folds,
                                                   length.out = length(unique_domains))),
                                 stringsAsFactors = FALSE)

  # Capture extra arguments for xgb
  dots <- list(...)

  if (verbose) {
    cat("Beginning cross-validation with population data\n")
  }

  # Set up parallel or sequential backend
  if (cpus > 1) {
    cl <- parallel::makeCluster(cpus)
    doSNOW::registerDoSNOW(cl)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    pb <- txtProgressBar(max = folds, style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    snow_opts <- list(progress = progress)
  } else {
    foreach::registerDoSEQ()
    if (verbose) {
      pb <- txtProgressBar(min = 0, max = folds, style = 3)
    }
    snow_opts <- list()
  }

  # Determine whether benchmark is being used
  use_benchmark <- !is.null(dots[["benchmark"]])

  if (cpus == 1) {
    # Sequential execution with a simple for loop for proper error reporting
    all_predictions_list <- vector("list", folds)

    for (fold in 1:folds) {
      # Identify held-out domains
      holdout_domains <- fold_assignment$Domain[fold_assignment$fold == fold]

      # Training sample: exclude held-out domains
      train_smp <- smp_data[!(smp_data[, domains] %in% holdout_domains), ]

      # Skip fold if no training data or no held-out domains
      if (nrow(train_smp) == 0 || length(holdout_domains) == 0) next

      # Build xgb arguments
      xgb_args <- list(
        fixed = fixed,
        smp_data = train_smp,
        smp_weights = smp_weights,
        pop_data = pop_data,
        pop_weights = pop_weights,
        domains = domains,
        sub_domains = sub_domains,
        transformation = transformation,
        bootstrap = FALSE,
        L = L,
        nrounds = nrounds,
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
        na.rm = na.rm,
        seed = seed,
        rescale_weights = rescale_weights,
        variance_y = variance_y,
        benchmark = NULL,
        benchmark_level = NULL,
        benchmark_weights = NULL
      )
      # Override with any user-supplied arguments (including benchmark args if provided)
      if (length(dots) > 0) {
        xgb_args[names(dots)] <- dots
      }

      fold_model <- suppressMessages(do.call(xgb, xgb_args))

      # Extract predictions for held-out domains
      fold_preds <- fold_model$ind
      holdout_preds <- fold_preds[as.character(fold_preds$Domain) %in% as.character(holdout_domains), ]

      # Use benchmarked estimates if benchmark is specified
      if ("Mean_bench" %in% colnames(holdout_preds)) {
        holdout_preds <- data.frame(Domain = holdout_preds$Domain,
                                    Predicted = holdout_preds$Mean_bench)
      } else {
        holdout_preds <- data.frame(Domain = holdout_preds$Domain,
                                    Predicted = holdout_preds$Mean)
      }

      all_predictions_list[[fold]] <- holdout_preds

      if (verbose) {
        setTxtProgressBar(pb, fold)
      }
    }

    all_predictions <- do.call(rbind, Filter(is.data.frame, all_predictions_list))

  } else {
    # Parallel execution with foreach
    all_predictions <- foreach::foreach(
      fold = 1:folds,
      .combine = function(...) {
        args <- list(...)
        args <- Filter(is.data.frame, args)
        if (length(args) == 0) return(NULL)
        do.call(rbind, args)
      },
      .packages = c("xgboost", "povmap"),
      .errorhandling = "remove",
      .options.snow = snow_opts
    ) %dopar% {

      # Identify held-out domains
      holdout_domains <- fold_assignment$Domain[fold_assignment$fold == fold]

      # Training sample: exclude held-out domains
      train_smp <- smp_data[!(smp_data[, domains] %in% holdout_domains), ]

      # Skip fold if no training data or no held-out domains
      if (nrow(train_smp) == 0 || length(holdout_domains) == 0) return(NULL)

      # Build xgb arguments
      xgb_args <- list(
        fixed = fixed,
        smp_data = train_smp,
        smp_weights = smp_weights,
        pop_data = pop_data,
        pop_weights = pop_weights,
        domains = domains,
        sub_domains = sub_domains,
        transformation = transformation,
        bootstrap = FALSE,
        L = L,
        nrounds = nrounds,
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
        na.rm = na.rm,
        seed = seed,
        rescale_weights = rescale_weights,
        variance_y = variance_y,
        benchmark = NULL,
        benchmark_level = NULL,
        benchmark_weights = NULL
      )
      # Override with any user-supplied arguments (including benchmark args if provided)
      if (length(dots) > 0) {
        xgb_args[names(dots)] <- dots
      }

      fold_model <- tryCatch({
        suppressMessages(do.call(xgb, xgb_args))
      }, error = function(e) {
        return(NULL)
      })

      if (is.null(fold_model)) return(NULL)

      # Extract predictions for held-out domains
      fold_preds <- fold_model$ind
      holdout_preds <- fold_preds[as.character(fold_preds$Domain) %in% as.character(holdout_domains), ]

      # Use benchmarked estimates if benchmark is specified
      if ("Mean_bench" %in% colnames(holdout_preds)) {
        holdout_preds <- data.frame(Domain = holdout_preds$Domain,
                                    Predicted = holdout_preds$Mean_bench)
      } else {
        holdout_preds <- data.frame(Domain = holdout_preds$Domain,
                                    Predicted = holdout_preds$Mean)
      }

      return(holdout_preds)
    }
  }

  if (verbose || cpus > 1) {
    close(pb)
    cat("\n")
  }

  # Check that we have results
  if (is.null(all_predictions) || 
      (is.data.frame(all_predictions) && nrow(all_predictions) == 0)) {
    stop("All cross-validation folds failed. Run with cpus=1 and verbose=TRUE to see errors.")
  }
  cv_predictions <- all_predictions
  colnames(cv_predictions) <- c("Domain", "Predicted")

  # Merge with direct estimates
  cv_results <- merge(cv_predictions, direct_estimates, by = "Domain")

  # Compute diagnostics on the back-transformed (original) scale.
  # All out-of-sample statistics are scored against the DIRECT estimate, which is
  # itself a noisy target, not against a known truth.
  mse_cv <- mean((cv_results$Direct - cv_results$Predicted)^2)
  mae_cv <- mean(abs(cv_results$Direct - cv_results$Predicted))
  # r2_cv: proportion of variance explained, exactly 1 - SSE/SST. Numerator and
  # denominator use matched denominators (both are sums). Penalizes bias and
  # scale compression and can be negative. Distinct from cor2_cv below.
  # (Previously this divided SSE by n and SST by n-1, understating the ratio by a
  # factor (n-1)/n and overstating r2_cv; corrected to matched sums.)
  sse_cv <- sum((cv_results$Direct - cv_results$Predicted)^2)
  sst_cv <- sum((cv_results$Direct - mean(cv_results$Direct))^2)
  r2_cv <- 1 - sse_cv / sst_cv
  # cor_cv: Pearson correlation between held-out prediction and direct estimate;
  # cor2_cv its square. Invariant to a linear rescaling of the prediction, so it
  # measures co-variation only and, unlike r2_cv, is not driven down by
  # compression of the prediction toward the mean.
  cor_cv <- cor(cv_results$Direct, cv_results$Predicted)
  cor2_cv <- cor_cv^2
  rank_cor_cv <- cor(cv_results$Direct, cv_results$Predicted, method = "spearman")

  # Skewness of domain-level outcomes
  skewness_domains <- mean(((cv_results$Direct - mean(cv_results$Direct)) /
                              sd(cv_results$Direct))^3)

  # Output
  result <- list(
    r2_cv = r2_cv,
    cor_cv = cor_cv,
    cor2_cv = cor2_cv,
    mse_cv = mse_cv,
    mae_cv = mae_cv,
    rank_cor_cv = rank_cor_cv,
    skewness_domains = skewness_domains,
    domain_results = cv_results
  )

  if (verbose) {
    cat("Cross-validation results (with population data):\n")
    cat("  (all scored against the noisy direct estimate, not truth)\n")
    cat("  R2 (1 - SSE/SST): ", round(r2_cv, 4), "\n")
    cat("  Correlation:      ", round(cor_cv, 4), "\n")
    cat("  Cor-squared:      ", round(cor2_cv, 4), "\n")
    cat("  MSE:              ", round(mse_cv, 6), "\n")
    cat("  MAE:              ", round(mae_cv, 4), "\n")
    cat("  Rank cor:         ", round(rank_cor_cv, 4), "\n")
    cat("  Skewness:         ", round(skewness_domains, 3), "\n")
  }

  # Generate predicted vs direct plot if filename is specified
  if (!is.null(plot_filename)) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
      warning("ggplot2 is required for plotting. Skipping plot.")
    } else {
      axis_max <- max(c(cv_results$Direct, cv_results$Predicted), na.rm = TRUE)
      axis_min <- min(c(cv_results$Direct, cv_results$Predicted), na.rm = TRUE)

      p <- ggplot2::ggplot(cv_results, ggplot2::aes(x = Direct, y = Predicted)) +
        ggplot2::geom_point(alpha = 0.4, size = 1.5) +
        ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
        ggplot2::coord_cartesian(xlim = c(axis_min, axis_max), ylim = c(axis_min, axis_max)) +
        ggplot2::labs(
          x = "Direct estimate (sample)",
          y = "Predicted (out-of-sample)",
          title = "Cross-validation: predicted vs. direct estimates",
          subtitle = paste0("R\u00B2 = ", round(r2_cv, 3),
                            ",  MAE = ", round(mae_cv, 4),
                            ",  Rank cor = ", round(rank_cor_cv, 3))
        ) +
        ggplot2::theme_minimal()

      ggplot2::ggsave(plot_filename, plot = p, width = 7, height = 6, dpi = 300)
      if (verbose) cat("  Plot saved to:", plot_filename, "\n")
      result$plot <- p
    }
  }

  class(result) <- "xgb_cv"
  return(result)
}
