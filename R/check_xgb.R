xgb_check1 <- function(transformation,
                       Y_smp,
                       X_smp,
                       X_pop,
                       smp_weights,
                       pop_weights,
                       conf_level,
                       domains,
                       sub_domains,
                       benchmark,
                       benchmark_level,
                       benchmark_weights,
                       benchmark_type){

  #if (!("dplyr" %in% .packages())) stop("dplyr package is required, please use library(dplyr) if installed")
  #require(dplyr)
  #require(xgboost)
  require(stats)

   if(!(transformation %in% c("no", "arcsin", "log","logistic","log.shift","ordernorm","poisson"))) stop("For transformation, please choose no, arcsin, logistic, log, log.shift, or poisson")

  if (transformation=="arcsin" | transformation=="logistic"){
    if(min(Y_smp)<0 | max(Y_smp)>1) stop("The outcome variable must be between 0 and 1 for arcsin or logistic transformations.")
  }

  if (transformation=="log"){
    if(min(Y_smp)<=0) stop("The outcome variable must be strictly greater than 0 for log transformations.")
  }

  #if(sum(is.na(Y_smp))>0) stop("There are missing values in the outcome variable.")

  #if(sum(is.na(X_smp))>0) stop("There are missing values in the independent variables in the sample dataset.")

  #if(sum(is.na(X_pop))>0) stop("There are missing values in the independent variables in the population dataset")

  if(sum(is.na(smp_weights))>0) stop("There are missing values in the sample weights.")

  if(sum(is.na(pop_weights))>0) stop("There are missing values in the population weights.")

  if(conf_level<=0 | conf_level>=1) stop("Please specify a confidence level between 0 and 1 (e.g. 0.95).")

  if (nrow(as.data.frame(Y_smp))!=nrow(X_smp)) stop("The number of rows in the outcome variable and independent variables are different.")

  if (length(which(colnames(X_smp)==paste0(domains)))==0) stop("The domain variable is not in the sample data.")

  if (length(which(colnames(X_pop)==paste0(domains)))==0) stop("The domain variable is not in the population data.")

  if (!is.character(domains)) stop("The domain name must be a character value.")

  if (!is.character(sub_domains)) stop("The subdomain name must be a character value.")

  #if (nrow(as.data.frame(smp_weights))!=nrow(X_smp)) stop("The number of rows in the sample weight column does not equal the number for the independent variables.")

  #if (is.null(smp_weights)==FALSE & length(smp_weights)!=nrow(X_smp)){
  #  stop("The length of weight variable does not equal the the number of rows in the sample data.")
  #}

  if (length(pop_weights)!=nrow(X_pop)) stop("The length of the population weight vector does not equal the number of rows in the population data.")

  #if (is.null(pop_weights)==FALSE & length(pop_weights)!=nrow(X_pop)){
  #  stop("Length of population weights must be the same as the number of rows of the population features.")
  #}
  if (!is.null(benchmark)) {
    if (!(is.numeric(benchmark) || is.character(benchmark) ||
          is.data.frame(benchmark))){
      stop(strwrap(prefix = " ", initial = "",
                   "For fixed value: Benchmark must be a named vector
                   containing the numeric benchmark value(s) and is of class
                   numeric. The names of the vector matchs to the chosen
                   indicators. \n For survey values: Benchmark must be a
                   vector of class character containing the names of the chosen
                   indicators."))
    }
    if (is.numeric(benchmark)) {
      if (!length(benchmark) %in% 1:2) {
        stop(strwrap(prefix = " ", initial = "",
                     "Benchmark must be a named vector containing the numeric
                     benchmark value(s) and is of class numeric. Benchmarking is
                     supplied for the Mean and the Head_Count ratio. Therefore,
                     the length of benchmark must be 1 or 2. The names of this
                     vector indicates whether the Mean, the Head_Count, or
                     both in which order are supplied."))
      }
      if (is.null(names(benchmark))) {
        stop(strwrap(prefix = " ", initial = "",
                     "Benchmark must be a named vector containing the numeric
                     benchmark value(s) and is of class numeric. Please provide
                     names."))
      }
      if (!length(benchmark) == length(names(benchmark))) {
        stop(strwrap(prefix = " ", initial = "",
                     "Benchmark must be a named vector containing the numeric
                     benchmark value(s) and is of class numeric. Each numeric must
                     be labeled. Therefore, benchmark and names(benchmark) have
                     the same length."))
      }
      if (!all(names(benchmark) %in% c("Mean", "Head_Count"))) {
        stop(strwrap(prefix = " ", initial = "",
                     "Benchmark must be a named vector containing the numeric
                     benchmark value(s) and is of class numeric. Benchmarking is
                     supplied for the Mean and the Head_Count ratio. Therefore,
                     the names must match with 'Mean' and 'Head_Count'."))
      }
      # if (!is.null(benchmark_weights)) {
      #   stop(strwrap(prefix = " ", initial = "",
      #                "For external benchmarking no benchmark weights can be
      #                used."))
      # }
    }
    if (is.character(benchmark)) {
      if(!length(benchmark) %in% 1:2) {
        stop(strwrap(prefix = " ", initial = "",
                     "Benchmark must be a vector of class character containing
                     the names of the chosen indicators. Benchmarking is
                     supplied for the Mean and the Head_Count ratio. Therefore,
                     the length of benchmark must be 1 or 2. The vector
                     indicates whether the Mean, the Head_Count, or both in
                     which order are supplied."))
      }
      if (!all(benchmark %in% c("Mean", "Head_Count"))) {
        stop(strwrap(prefix = " ", initial = "",
                     "Benchmark must be a vector of class character containing
                     the names of the chosen indicators. Benchmarking is
                     supplied for the Mean and the Head_Count ratio. Therefore,
                     it must match with 'Mean' and 'Head_Count'."))
      }
      if (is.null(weights)) {
        stop(strwrap(prefix = " ", initial = "",
                     "The argument benchmark indicates that it is benchmarked
                     with the survey data. Please provide weights through the
                     argument weights."))
      }
    }
    if (is.data.frame(benchmark)) {
      if (is.null(benchmark_level)) {
        stop(strwrap(prefix = " ", initial = "",
                     "As the input in benchmark is a data.frame. Fixed benchmark
                     values are used at a lower level. Please give the name
                     of this variable in the sample and population data
                     by the argument benchmark_level."))
      }
      if (!length(benchmark) %in% 2:3) {
        stop(strwrap(prefix = " ", initial = "",
                     "Benchmark must be a data.frame composed of a variable
                     of class character containing the domain names at which the
                     benchmarkaing is performed and variable(s) with
                     benchmark value(s) of class numeric. Benchmarking is
                     supplied for the Mean and the Head_Count ratio. Therefore,
                     the names of the data.frame must match for the first
                     variable the benchmark_level and for the other(s) to Mean
                     and Head_Count."))
      }
      if (is.null(names(benchmark))) {
        stop(strwrap(prefix = " ", initial = "",
                     "Benchmark must be a data.frame composed of a variable
                     of class character containing the domain names at which the
                     benchmarkaing is performed and variable(s) with
                     benchmark value(s) of class numeric. Benchmarking is
                     supplied for the Mean and the Head_Count ratio. Therefore,
                     the names of the data.frame must match for the first
                     variable the benchmark_level and for the other(s) to Mean
                     and Head_Count. Please provide names."))
      }
      if (!length(benchmark) == length(names(benchmark))) {
        stop(strwrap(prefix = " ", initial = "",
                     "Benchmark must be a data.frame composed of a variable
                     of class character containing the domain names at which the
                     benchmarkaing is performed and variable(s) with
                     benchmark value(s) of class numeric. Benchmarking is
                     supplied for the Mean and the Head_Count ratio. Therefore,
                     the names of the data.frame must match for the first
                     variable the benchmark_level and for the other(s) to Mean
                     and Head_Count. Each variable in the data.frame must
                     be labeled."))
      }
      if (!all(names(benchmark)[-1] %in% c("Mean", "Head_Count"))) {
        stop(strwrap(prefix = " ", initial = "",
                     "Benchmark must be a data.frame composed of a variable
                     of class character containing the domain names at which the
                     benchmarkaing is performed and variable(s) with
                     benchmark value(s) of class numeric. Benchmarking is
                     supplied for the Mean and the Head_Count ratio. Therefore,
                     the names of the data.frame must match for the first
                     variable the benchmark_level and for the other(s) to Mean
                     and Head_Count. No other names are possible."))
      }
      if (names(benchmark)[1] != benchmark_level) {
        stop(strwrap(prefix = " ", initial = "",
                     "Benchmark must be a data.frame composed of a variable
                     of class character containing the domain names at which the
                     benchmarkaing is performed and variable(s) with
                     benchmark value(s) of class numeric. Benchmarking is
                     supplied for the Mean and the Head_Count ratio. Therefore,
                     the names of the data.frame must match for the first
                     variable the benchmark_level and for the other(s) to Mean
                     and Head_Count. The name of the first variable indicataing
                     the domains of the benchmark_level does not match to the
                     argument benchmark_level."))
      }
      # if (!is.null(benchmark_weights)) {
      #   stop(strwrap(prefix = " ", initial = "",
      #                "For external benchmarking no benchmark weights can be
      #                used."))
      # }
    }
  }

  if (!is.null(benchmark_level)) {
    if (!benchmark_level %in% colnames(X_pop)) {
      stop(strwrap(prefix = " ", initial = "",
                   paste0("The variable ",benchmark_level, " specified as the benchmark_level is not contained in the population data")))
    }
  }


  if (benchmark_type != "ratio" && benchmark_type != "raking" && benchmark_type != "ratio_complement" && benchmark_type != "ratio_bound" && benchmark_type != "logit_raking") {
    stop(strwrap(prefix = " ", initial = "",
                 "The benchmark version of ebp is only available with
                   'raking', 'ratio', 'ratio_complement', 'ratio_bound', and 'logit_raking'."))
  }

  if (benchmark_type == "ratio_complement" && is.data.frame(benchmark))  {
    if (max(benchmark[["Head_Count"]])>1 | max(benchmark[["Mean"]])>1) {
      stop(strwrap(prefix = " ", initial = "",
                   "When benchmarking with ratio_complement, the target values must lie between 0 and 1."))
    }
  }


  if (is.null(benchmark) && benchmark_type != "ratio") {
    stop(strwrap(prefix = " ", initial = "",
                 "A benchmark type is provided, but no benchmark value.
                   Please provide the argument 'benchmark' within the
                   function."))

}







}

xgb_check2 <- function(transformation,
                       Y_smp,
                       X_smp,
                       smp_weights,
                       domains,
                       cluster){




 if(!(transformation %in% c("no", "arcsin", "log","logistic","log.shift","ordernorm","poisson"))) stop("For transformation, please choose no, arcsin, logistic, log, log.shift, or poisson")

  if (transformation=="arcsin"){
    if(min(Y_smp)<0 | max(Y_smp)>1) stop("The outcome variable must be between 0 and 1 for arcsin transformations.")
  }

  if (transformation=="log"){
    if(min(Y_smp)<=0) stop("The outcome variable must be strictly greater than 0 for log transformations.")
  }

  if (transformation=="poisson"){
    if(min(Y_smp)<0) stop("The outcome variable must be non-negative for the poisson transformation.")
  }

  if(sum(is.na(Y_smp))>0) stop("There are missing values in the outcome variable.")

  if(sum(is.na(X_smp[,domains]))>0) stop("There are missing values in the domain variable in the sample dataset.")

  if(sum(is.na(smp_weights))>0) stop("There are missing values in the sample weights.")

  if (length(colnames(Y_smp))>1) stop("The outcome variable must be a vector or have just one column.")

  if (nrow(as.data.frame(Y_smp))!=nrow(X_smp)) stop("The lengths of the outcome variable and independent variables are different.")

  if (length(which(colnames(X_smp)==paste0(cluster)))==0 & cluster!=domains) stop("The cluster variable is not in the sample data.")

  if (length(which(colnames(X_smp)==paste0(domains)))==0) stop("The domain variable is not in the sample data.")

  if (nrow(as.data.frame(smp_weights))!=nrow(as.data.frame(Y_smp))) stop("The length of the weight variable does not equal the number of observations in the outcome variable.")

  if (is.null(as.data.frame(smp_weights))==FALSE & nrow(as.data.frame(smp_weights))!=nrow(X_smp)){
    stop("The length of weight variable does not equal the the number of rows in the sample data.")
  }

  if (!is.character(domains)) stop("The domain name must be a character value.")

  if (!is.character(cluster)) stop("The cluster name must be a character value.")

}

#' Reject arguments that fell into `...` without being used
#'
#' `xgb_tune()` forwards `...` to `xgboost::xgb.train()`, so an argument that is
#' not a formal of `xgb_tune()` and not an xgboost parameter is silently
#' discarded. That is how a caller written against a newer signature can run
#' against an older library and get a result with no error: when `search`,
#' `n_iter`, `bounds` and `seed` were added, a call supplying them to a build
#' that predated them evaluated the 32-point default grid -- every point at
#' `eta = 0.3`, `nround` in {150, 300} -- and returned it looking like a random
#' search, with only a suspiciously short runtime to give it away.
#'
#' This check makes that failure loud. It cannot repair older builds, which have
#' no check in them; it protects every version from here on, which is the case
#' where a caller is upgraded before the library.
#'
#' @param dots the result of `list(...)` in the calling function.
#' @param fn_formals names of the calling function's formals, used to spot a
#'   near-miss spelling of a real argument (`fold` for `folds`).
#' @param allow character vector of additional names to accept.
#' @param strict if `TRUE` (default) an unrecognised name is an error; if
#'   `FALSE` it is a warning. The escape hatch exists so that a genuinely new
#'   xgboost parameter, added upstream after this list was written, can never
#'   block a legitimate call.
#' @keywords internal
#' @noRd
.check_xgb_dots <- function(dots, fn_formals = character(0),
                            allow = character(0), strict = TRUE) {
  if (!length(dots)) return(invisible(TRUE))

  nms <- names(dots)
  if (is.null(nms) || any(!nzchar(nms))) {
    stop("xgb_tune(): ", sum(is.null(nms) | !nzchar(nms)),
         " argument(s) were passed to `...` without a name. `...` is forwarded ",
         "to xgboost::xgb.train(), which needs names, so an unnamed value is ",
         "discarded. Name them, or remove them.", call. = FALSE)
  }

  ## Documented xgboost training parameters plus the xgb.train arguments a
  ## caller may reasonably want to reach. Deliberately generous: the point is to
  ## catch a name that belongs to no API at all, not to police xgboost.
  known <- c(
    ## general
    "booster", "device", "nthread", "verbosity", "validate_parameters", "seed_per_iteration",
    ## tree booster
    "base_score", "objective", "eval_metric", "tree_method", "grow_policy",
    "max_leaves", "max_bin", "sampling_method", "monotone_constraints",
    "interaction_constraints", "num_parallel_tree", "scale_pos_weight",
    "sketch_eps", "refresh_leaf", "process_type", "updater", "top_k",
    "colsample_bylevel", "colsample_bynode", "colsample_bytree",
    "max_cat_to_onehot", "max_cat_threshold", "multi_strategy",
    ## objective-specific
    "tweedie_variance_power", "huber_slope", "quantile_alpha",
    "aft_loss_distribution", "aft_loss_distribution_scale", "num_class",
    "disable_default_eval_metric",
    ## xgb.train arguments
    "watchlist", "evals", "obj", "feval", "custom_metric", "maximize",
    "print_every_n", "callbacks", "xgb_model", "save_period", "save_name",
    "missing", "weight", "params",
    allow)

  unknown <- setdiff(nms, known)
  if (!length(unknown)) return(invisible(TRUE))

  ## A name that nearly matches a real formal is almost always a typo or a
  ## signature change, so say which one it looks like rather than only that it
  ## was not recognised.
  hint <- character(0)
  if (length(fn_formals)) {
    fn_formals <- setdiff(fn_formals, "...")
    for (u in unknown) {
      d <- utils::adist(u, fn_formals, ignore.case = TRUE)[1, ]
      near <- fn_formals[d <= max(1L, floor(nchar(u) / 4))]
      if (length(near)) hint <- c(hint, sprintf("  '%s' -- did you mean %s?",
                                                u, paste(sQuote(near), collapse = " or ")))
    }
  }

  msg <- paste0(
    "xgb_tune(): unrecognised argument(s) in `...`: ",
    paste(sQuote(unknown), collapse = ", "), ".\n",
    "`...` is forwarded to xgboost::xgb.train(), so these would be silently ",
    "discarded rather than used.\n",
    if (length(hint)) paste0(paste(hint, collapse = "\n"), "\n") else "",
    "If this is an argument of a NEWER xgb_tune() than the one installed ",
    "(installed version ", utils::packageVersion("povmap"),
    "), the call would otherwise have run with defaults and returned a result ",
    "that looks correct.\n",
    "If it is a legitimate xgboost parameter this list does not know about, ",
    "pass it inside `params = list(...)`.")
  if (isTRUE(strict)) stop(msg, call. = FALSE) else warning(msg, call. = FALSE)
  invisible(FALSE)
}
