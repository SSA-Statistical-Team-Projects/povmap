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