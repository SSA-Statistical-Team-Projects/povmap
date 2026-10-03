mse_povmap <- function(object, indicator = "all", CV = FALSE) {
  if (inherits(object, "ell") | (inherits(object, "xgb"))) {
     object$MSE <- object$var
   }
  
   if (is.null(object$MSE) && CV == TRUE) {
    stop(strwrap(prefix = " ", initial = "",
                 "No MSE estimates in povmap object: arguments MSE and CV have to
                 be FALSE or a new povmap object with variance/MSE needs to be
                 generated."))
  }
  if ((ncol(object$ind) == 11) && any(indicator == "Custom" |
    indicator == "custom")) {
    stop(strwrap(prefix = " ", initial = "",
                 "No individual indicators are defined. Either select other
                 indicators or define custom indicators and generate a new povmap
                 object. See also help(ebp)."))
  }

  # Calculation of CVs
  if (inherits(object, "fh")) {
    object$MSE <- object$MSE[, c("Domain", "Direct", "FH")]
    object$ind <- object$ind[, c("Domain", "Direct", "FH")]
  }
  object <- align_precision(object)
  all_cv <- sqrt(object$MSE[, -1, drop = FALSE]) / object$ind[, -1, drop = FALSE]

  if (any(indicator == "Quantiles") || any(indicator == "quantiles")) {
    indicator <- c(
      indicator[!(indicator == "Quantiles" ||
        indicator == "quantiles")],
      "Quantile_10", "Quantile_25", "Median",
      "Quantile_75", "Quantile_90"
    )
  }
  if (any(indicator == "poverty") || any(indicator == "Poverty")) {
    indicator <- c(
      indicator[!(indicator == "poverty" ||
        indicator == "Poverty")],
      "Head_Count", "Poverty_Gap"
    )
  }
  if (any(indicator == "inequality") || any(indicator == "Inequality")) {
    indicator <- c(
      indicator[!(indicator == "inequality" ||
        indicator == "Inequality")],
      "Gini", "Quintile_Share"
    )
  }
  if (any(indicator == "custom") || any(indicator == "Custom")) {
    indicator <- c(
      indicator[!(indicator == "custom" | indicator == "Custom")],
      colnames(object$ind[-c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11)])
    )
  }

  if (any(indicator == "all") || any(indicator == "All")) {
    ind <- object$MSE
    ind_cv <- data.frame(Domain = object$MSE[, 1], all_cv)
    ind_name <- "All indicators"
  } else if (any(indicator == "fh") || any(indicator == "FH")) {
    ind <- object$MSE[, c("Domain", "FH")]
    ind_cv <- data.frame(Domain = object$MSE[, 1], all_cv)
    ind_name <- "Fay-Herriot estimates"
  } else if (any(indicator == "Direct") || any(indicator == "direct")) {
    ind <- object$MSE[, c("Domain", "Direct")]
    ind_cv <- data.frame(Domain = object$MSE[, 1], all_cv)
    ind_name <- "Direct estimates used in Fay-Herriot approach"
  } else {
    selection <- colnames(object$MSE[-1]) %in% indicator
    ind <- object$MSE[, c(TRUE, selection)]
    ind_cv <- data.frame(Domain = object$MSE[, 1], all_cv[, selection])
    colnames(ind_cv) <- colnames(ind)
    ind_name <- paste(unique(indicator), collapse = ", ")
  }

  if (CV == FALSE) {
    mse_povmap <- list(ind = ind, ind_name = ind_name)
  } else {
    mse_povmap <- list(ind = ind, ind_cv = ind_cv, ind_name = ind_name)
  }

  class(mse_povmap) <- "mse.povmap"

  return(mse_povmap)
}

# Pair each precision (MSE or variance) column with the point estimate it measures, by
# NAME rather than by position, before CV = sqrt(precision) / point is computed. Used by
# mse_povmap() and mse_emdi().
#   * xgb (and megb) store the variance of the benchmarked estimate in $var as
#     "Var_bench"; it measures $ind$Mean_bench, so it is renamed "Mean_bench" here.
#     $var$Mean measures $ind$Mean.
#   * A point column with no precision counterpart gets no CV, instead of being divided
#     by whichever precision column sits in the same position. Since c552f02, xgb's $ind
#     also carries Mean_agg (the aggregate-then-back-transform point, kept only to
#     reconcile with earlier versions). No bootstrap variance is computed for it, so it
#     gets no _Var or _CV column.
#   * Rows are matched on Domain.
# A precision column with no point counterpart means the object is inconsistent: stop.
align_precision <- function(object) {
  prec <- object$MSE
  pt <- object$ind
  names(prec)[names(prec) == "Var_bench"] <- "Mean_bench"
  orphan <- setdiff(names(prec)[-1], names(pt)[-1])
  if (length(orphan) > 0) {
    stop("No point estimate matches the precision column(s) ",
         paste(orphan, collapse = ", "), "; cannot compute CVs.")
  }
  rows <- match(as.character(prec[[1]]), as.character(pt[[1]]))
  if (anyNA(rows)) {
    stop("Some domains in the precision estimates have no point estimate.")
  }
  object$MSE <- prec
  object$ind <- pt[rows, c(names(pt)[1], names(prec)[-1]), drop = FALSE]
  rownames(object$ind) <- NULL
  object
}
