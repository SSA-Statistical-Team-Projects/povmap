# Internal documentation -------------------------------------------------------

# Benchmark function for the EBP with a national value--------------------------

# This function is called within the EBP-function and the agruments benchmark
# and benchmark_type are documented there.



# Benchmark function for the EBP with a variable domain value ------------------

# This function is called within the EBP-function and the agruments benchmark
# and benchmark_type are documented there.

benchmark_xgb_level <- function (point_estim, framework, fixed, benchmark,
                                 benchmark_type, benchmark_level,benchmark_weights = NULL) {



  if (!is.numeric(benchmark)) {
    benchmark_ <- collapse::fmean(framework$smp_data[[paste0(fixed[2])]],g=framework$smp_data[,benchmark_level],w=framework$smp_data[,framework$smp_weights])
    benchmark_ <- data.frame(benchmark_,names(benchmark_))
    colnames(benchmark_) <- c("Mean",benchmark_level)
  }
    else {
  benchmark_ <- as.data.frame(benchmark)
  benchmark_[,benchmark_level] <- names(benchmark)
  colnames(benchmark_) <- c("Mean",benchmark_level)
    }


  benchmark_df <- data.frame(point_estim$ind$Mean, unique(framework$pop_data[[framework$domains]]),1:length(point_estim$ind$Mean))
  colnames(benchmark_df) <- c("point_estim",framework$domains,"order")
  crosswalk <- data.frame(unique(framework$pop_data[,c(framework$domains,benchmark_level)]))


  colnames(crosswalk) <- c(framework$domains,benchmark_level)
  benchmark_df <- merge(benchmark_df,crosswalk,by=framework$domains,all.x=T,sort=F)
  benchmark_df <- merge(benchmark_df,benchmark_,by=benchmark_level,all.x=T,sort=F)
  # unfortunately because it is a many to one merge the sort order is not preserved, so we need to reorder by hand
  benchmark_df <- benchmark_df[order(benchmark_df$order),]
  popwts <- data.frame(collapse::fsum(framework$pop_data[,framework$pop_weights],g=framework$pop_data[,framework$domains]))
  population_lga2 <- collapse:::fsum(framework$pop_data$population,g=framework$pop_data$lgacode)
  colnames(popwts) <- framework$pop_weights
  benchmark_df <- data.frame(benchmark_df,popwts)
  weighted_pe <- data.frame(collapse::fmean(benchmark_df$point_estim,g=benchmark_df[,benchmark_level],w=benchmark_df[,framework$pop_weights],TRA=1))
  colnames(weighted_pe) <- "weighted_pe"
  benchmark_df <- data.frame(benchmark_df,weighted_pe)

  if (benchmark_type == "ratio" | benchmark_type=="ratio_bound") {
      #ratio benchmarking
        benchmark_df$xgb_bench <- benchmark_df$point_estim*(benchmark_df$Mean/benchmark_df$weighted_pe)
        # For benchmark levels with no survey data, do no benchmarking
        benchmark_df$xgb_bench[is.na(benchmark_df$Mean)]<-benchmark_df$point_estim[is.na(benchmark_df$Mean)]
  }
  if (benchmark_type=="ratio_complement") {
          factors_complement <- (1-benchmark_df$Mean)/(1-benchmark_df$weighted_pe)
          benchmark_df$xgb_bench <- 1-(1-benchmark_df$point_estim)*factors_complement
          benchmark_df$xgb_bench[is.na(benchmark_df$Mean)]<-benchmark_df$point_estim[is.na(benchmark_df$Mean)]
  }
  if (benchmark_type=="logit_raking") {
    eps <- 1e-7
    pe_clamped <- pmax(eps, pmin(1 - eps, benchmark_df$point_estim))
    logit_pe <- log(pe_clamped / (1 - pe_clamped))
    benchmark_df$xgb_bench <- benchmark_df$point_estim  # initialize
    levels <- unique(benchmark_df[, benchmark_level])
    for (lev in levels) {
      idx <- benchmark_df[, benchmark_level] == lev
      target <- benchmark_df$Mean[idx][1]
      if (is.na(target)) {
        benchmark_df$xgb_bench[idx] <- benchmark_df$point_estim[idx]
        next
      }
      # Boundary targets: logit raking cannot reach exactly 0 or 1
      if (target <= eps) {
        benchmark_df$xgb_bench[idx] <- 0
        next
      }
      if (target >= 1 - eps) {
        benchmark_df$xgb_bench[idx] <- 1
        next
      }
      wts <- benchmark_df[idx, framework$pop_weights]
      logit_i <- logit_pe[idx]
      # Find additive shift c such that weighted mean of logit_inv(logit_i + c) = target
      obj <- function(c) {
        adjusted <- exp(logit_i + c) / (1 + exp(logit_i + c))
        weighted.mean(adjusted, w = wts) - target
      }
      # Check if shift is needed
      if (abs(obj(0)) < 1e-12) {
        benchmark_df$xgb_bench[idx] <- pe_clamped[idx]
        next
      }
      sol <- uniroot(obj, interval = c(-50, 50), tol = 1e-10)
      adjusted <- exp(logit_i + sol$root) / (1 + exp(logit_i + sol$root))
      benchmark_df$xgb_bench[idx] <- adjusted
    }
  }
  if (benchmark_type=="ratio_bound") {
    max <- collapse::fmax(benchmark_df$xgb_bench,g=benchmark_df[,benchmark_level],TRA=1)
    factors_complement <- (1-benchmark_df$Mean)/(1-benchmark_df$weighted_pe)
    benchmark_df$xgb_bench[max>1] <- 1-(1-benchmark_df$point_estim[max>1])*factors_complement[max>1]
    benchmark_df$xgb_bench[is.na(benchmark_df$Mean)]<-benchmark_df$point_estim[is.na(benchmark_df$Mean)]
  }

# check
#a <- collapse:::fmean(benchmark_df$xgb_bench,g=benchmark_df$state_,w=benchmark_df$population)

  if (is.list(point_estim)) {
    point_estim_bench <- data.frame(benchmark_df$point_estim, benchmark_df$xgb_bench)
  } else {
    point_estim_bench <- as.matrix(data.frame(benchmark_df$point_estim, benchmark_df$xgb_bench))
  }
  colnames(point_estim_bench)[1:2] <- c("Mean","Mean_bench")
  return(point_estim_bench)
}
