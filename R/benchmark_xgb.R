# Internal documentation -------------------------------------------------------

# Benchmark function for the EBP with a national value--------------------------

# This function is called within the EBP-function and the agruments benchmark
# and benchmark_type are documented there.



# Benchmark function for the EBP with a variable domain value ------------------

# This function is called within the EBP-function and the agruments benchmark
# and benchmark_type are documented there.

benchmark_xgb_level <- function (point_estim, framework, fixed, benchmark,
                                 benchmark_type, benchmark_level) {

  
  if (!is.numeric(benchmark)) {
    superarea <- benchmark_level 
    benchmark_ <- collapse::fmean(framework$smp_data[[paste0(fixed[2])]],g=framework$smp_data[,superarea],w=framework$smp_data[,framework$smp_weights])
    benchmark_ <- data.frame(benchmark_,names(benchmark_))
    colnames(benchmark_) <- c("Mean",superarea)
  }
    else {
  superarea <- names(benchmark)[1]
  benchmark_ <- benchmark 
    }
  benchmark_df <- data.frame(point_estim$ind$Mean, unique(framework$pop_data[[framework$domains]]),1:length(point_estim$ind$Mean)) 
  colnames(benchmark_df) <- c("point_estim",framework$domains,"order")
  crosswalk <- data.frame(unique(framework$pop_data[,c(framework$domains,superarea)]))
  
  
  colnames(crosswalk) <- c(framework$domains,superarea)
  benchmark_df <- merge(benchmark_df,crosswalk,by=framework$domains,all.x=T,sort=F)
  benchmark_df <- merge(benchmark_df,benchmark_,by=superarea,all.x=T,sort=F)
  # unfortunately because it is a many to one merge the sort order is not preserved, so we need to reorder by hand 
  benchmark_df <- benchmark_df[order(benchmark_df$order),]
  popwts <- data.frame(collapse::fsum(framework$pop_data[,framework$pop_weights],g=framework$pop_data[,framework$domains]))
  colnames(popwts) <- framework$pop_weights
  benchmark_df <- data.frame(benchmark_df,popwts)
  weighted_pe <- data.frame(collapse::fmean(benchmark_df$point_estim,g=benchmark_df[,superarea],w=benchmark_df[,framework$pop_weights],TRA=1))
  colnames(weighted_pe) <- "weighted_pe"                         
  benchmark_df <- data.frame(benchmark_df,weighted_pe)
                              
  if (benchmark_type == "ratio" | benchmark_type=="ratio_bound") {
      #ratio benchmarking 
        benchmark_df$xgb_bench <- benchmark_df$point_estim*(benchmark_df$Mean/benchmark_df$weighted_pe)
  }
  if (benchmark_type=="ratio_complement") {
          factors_complement <- (1-benchmark_df$Mean)/(1-benchmark_df$weighted_pe)
          benchmark_df$xgb_bench <- 1-(1-benchmark_df$point_estim)*factors_complement
  }
  if (benchmark_type=="ratio_bound") {
    max <- collapse::fmax(benchmark_df$xgb_bench,g=benchmark_df[,superarea],TRA=1)
    factors_complement <- (1-benchmark_df$Mean)/(1-benchmark_df$weighted_pe)
    benchmark_df$xgb_bench[max>1] <- 1-(1-benchmark_df$point_estim[max>1])*factors_complement[max>1]  
  }  
   
 

  if (is.list(point_estim)) {
    point_estim_bench <- data.frame(benchmark_df$point_estim, benchmark_df$xgb_bench)
  } else {
    point_estim_bench <- as.matrix(data.frame(benchmark_df$point_estim, benchmark_df$xgb_bench))
  }
  colnames(point_estim_bench)[1:2] <- c("Mean","Mean_bench")
  return(point_estim_bench)
}
