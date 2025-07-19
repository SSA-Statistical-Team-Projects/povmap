xgb_check1 <- function(transformation,
                       Y_smp,
                       X_smp,
                       X_pop,
                       smp_weights,
                       pop_weights,
                       conf_level,
                       domains,
                       sub_domains){

  if(!(transformation %in% c("no", "arcsin", "log"))) stop("For transformation, please choose no, arcsin, or log.")

  if (transformation=="arcsin"){
    if(min(Y_smp)<0 | max(Y_smp)>1) stop("The outcome variable must be between 0 and 1 for arcsin transformations.")
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

  if (length(Y_smp)!=nrow(X_smp)) stop("The lengths of the outcome variable and independent variables are different.")

  if (length(which(colnames(X_smp)==paste0(domains)))==0) stop("The domain variable is not in the sample data.")

  if (length(which(colnames(X_pop)==paste0(domains)))==0) stop("The domain variable is not in the population data.")

  if (!is.character(domains)) stop("The domain name must be a character value.")

  if (!is.character(sub_domains)) stop("The subdomain name must be a character value.")

  if (length(smp_weights)!=nrow(X_smp)) stop("The length of the sample weight vector does not equal the number of independent variables.")

  #if (is.null(smp_weights)==FALSE & length(smp_weights)!=nrow(X_smp)){
  #  stop("The length of weight variable does not equal the the number of rows in the sample data.")
  #}

  if (length(pop_weights)!=nrow(X_pop)) stop("The length of the population weight vector does not equal the number of independent variables.")

  #if (is.null(pop_weights)==FALSE & length(pop_weights)!=nrow(X_pop)){
  #  stop("Length of population weights must be the same as the number of rows of the population features.")
  #}
}

xgb_check2 <- function(transformation,
                       Y_smp,
                       X_smp,
                       smp_weights,
                       domains,
                       cluster){

  if(!(transformation %in% c("no", "arcsin", "log"))) stop("For transformation, please choose no, arcsin, or log.")

  if (transformation=="arcsin"){
    if(min(Y_smp)<0 | max(Y_smp)>1) stop("The outcome variable must be between 0 and 1 for arcsin transformations.")
  }

  if (transformation=="log"){
    if(min(Y_smp)<=0) stop("The outcome variable must be strictly greater than 0 for log transformations.")
  }

  if(sum(is.na(Y_smp))>0) stop("There are missing values in the outcome variable.")

  if(sum(is.na(X_smp))>0) stop("There are missing values in the independent variables in the sample dataset.")

  if(sum(is.na(smp_weights))>0) stop("There are missing values in the sample weights.")

  if (length(colnames(Y_smp))>1) stop("The outcome variable must be a vector or have just one column.")

  if (nrow(Y_smp)!=nrow(X_smp)) stop("The lengths of the outcome variable and independent variables are different.")

  if (length(which(colnames(X_smp)==paste0(cluster)))==0) stop("The cluster variable is not in the sample data.")

  if (length(which(colnames(X_smp)==paste0(domains)))==0) stop("The domain variable is not in the sample data.")

  if (length(smp_weights)!=nrow(Y_smp)) stop("The length of the weight variable does not equal the number of observations in the outcome variable.")

  if (is.null(smp_weights)==FALSE & length(smp_weights)!=nrow(X_smp)){
    stop("The length of weight variable does not equal the the number of rows in the sample data.")
  }

  if (!is.character(domains)) stop("The domain name must be a character value.")

  if (!is.character(cluster)) stop("The cluster name must be a character value.")

}
