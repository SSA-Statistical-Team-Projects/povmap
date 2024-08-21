point_estim_hdp <- function(framework,
                           fixed,
                           transformation,
                           interval,
                           L,
                           keep_data = FALSE,
                           Ydump = NULL){
  # Transformation of Data------------------------------------------------------
  # Optimal parameter function returns the minimum of the optimization
  # functions from generic_opt; the minimum is the optimal lambda.
  # The function can be found in the script optimal_parameter.R
  optimal_lambda <- optimal_parameter(
    generic_opt = generic_opt,
    fixed = fixed,
    smp_data = framework$smp_data,
    smp_domains = framework$smp_domains,
    transformation = transformation,
    interval = interval,
    framework = framework
  )
  
  # Data_transformation function returns transformed data and shift parameter.
  # The function can be found in the script transformation_functions.R
  transformation_par <- data_transformation(
    fixed = fixed,
    smp_data = framework$smp_data,
    transformation = transformation,
    lambda = optimal_lambda
  )
  shift_par <- transformation_par$shift
  
  # Model estimation, model parameter and parameter of generating model---------
  model.f <- model.frame(fixed, data = transformation_par$transformed_data)
  ys <- as.numeric(model.response(model.f))
  xs <- model.matrix(fixed, data = transformation_par$transformed_data)[ , , drop = FALSE]
  area.s <- transformation_par$transformed_data[[framework$smp_domains]]
  p <- ncol(xs)
  area_id <- sort(unique(area.s))
  #step 1: MQ function to get the q_ij
  if(is.null(framework$weights)){
    var.weights = rep(1, nrow(xs))
  }
  else{var.weights = framework$weights}
  mod <- QRLM(y = ys, x = xs, q = framework$Q, var.weights = var.weights, maxit = framework$maxit, acc = framework$tol, k = framework$k_b)
  qo <- matrix(c(gridfitinter(ys, mod$fitted.values, mod$q.values)), nrow = framework$N_smp, ncol=1)
  qmat <- matrix(c(qo, area.s), nrow = framework$N_smp, ncol = 2)
  #step 2: LBP to get the Q_i
  if(framework$method == "LBP"){
  tmp.Xmean <- rep(1,framework$N_dom_smp)
  tmp.Popn <- framework$n_pop[framework$dist_obs_dom]
  Xmean <- data.frame(area_id, tmp.Xmean)
  Popn  <- data.frame(area_id, tmp.Popn)
  tmp.dta<-data.frame(qmat[,1],area.s)
  names(tmp.dta) <- c("qij", "area.s")
  LBP <- eblupBHF(formula = qij~1, dom = area.s, meanxpop = Xmean, popnsize = Popn,
                method = "REML", data=tmp.dta)
  Qi <- LBP$eblup$eblup
  Qi <- ifelse(Qi>1, 0.999, Qi)
  Qi <- ifelse(Qi<0, 0.001, Qi)
  }
  if(framework$method == "Naive"){Qi <- tapply(qmat[,1], qmat[,2], mean)}
  #step 3: estimate Beta_Q_i, sigma^2_v_i and sigma^2_e_i for in-sample domains
  # mixed_model <- mqre_g(qtl=Qi,y=ys,x=xs,group=area.s,tol=framework$tol,maxit=framework$maxit,k=framework$k_b,
               # k_sigma_u=framework$k_sigma_u,k_sigma_e=framework$k_sigma_e)
  mixed_model <- mqre_g(fixed = fixed,
                        framework = framework,
                        transformation_par = transformation_par,
                        Qi = Qi,
                        var.weights = var.weights,
                        maxit = framework$maxit, 
                        acc = framework$tol,
                        k_b = framework$k_b, 
                        k_sigma_e = framework$k_sigma_e,
                        k_sigma_u = framework$k_sigma_u)

  #regression coefficients for out-of-sample area via using tau = 0.5
  mod.out <- QRLM(y = ys,
                  x = xs,
                  q = 0.5,
                  var.weights = var.weights,
                  maxit = framework$maxit,
                  acc = framework$tol,
                  k = framework$k_b)
  betas.out <- mod.out$coef
  # ends parameter estimation
  
  # Function model_par extracts the needed parameters  
  est_par <- model_par_hdp(
    mixed_model = mixed_model,
    transformation_par = transformation_par,
    framework = framework,
    fixed = fixed,
    betas.out = betas.out
  )
  
  # Function gen_model calculates the parameters in the generating model.
  gen_par <- gen_model_hdp(
    model_par = est_par,
    fixed = fixed,
    framework = framework
  )
  
  # Monte-Carlo approximation --------------------------------------------------
  if (inherits(framework$threshold, "function")) {
    framework$threshold <-
      framework$threshold(
        y = as.numeric(framework$smp_data[[paste0(fixed[2])]])
      )
  }
  
  # The monte-carlo function returns a data frame of desired indicators.
  indicator_prediction <- monte_carlo(
    L = L,
    framework = framework,
    model_par = est_par,
    gen_model = gen_par,
    fixed = fixed,
    transformation = transformation,
    lambda = optimal_lambda,
    shift = shift_par
  )
  
  return(list(
    ind = indicator_prediction,
    optimal_lambda = optimal_lambda,
    shift_par = shift_par,
    model_par = est_par,
    gen_model = gen_par,
    model = mixed_model
  ))
}

#Functions for internal use-----------------------------------------------------

#Finding the Quantile Orders by Linear Interpolation----------------------------
# Assumes that "zerovalinter" function has been already loaded
"gridfitinter" <- function(y,expectile,Q)
  # computing of the expectile-order of each observation of y by interpolation
{
  nq <- length(Q)
  diff <- y %*% t(as.matrix(rep(1, nq))) - expectile
  vectordest <- apply(diff, 1, zerovalinter,Q)
  #print(vectordest)
  #qord<-list(ord=c(vectordest))
  #qord
}
# COMPUTING OF THE QUANTILE-ORDERS--------------------------------------------
"zerovalinter" <- function(y, x)
{
  if(min(y) > 0) {
    xmin <- x[y == min(y)]
    if(length(xmin) > 0)
      xmin <- xmin[length(xmin)]
    xzero <- xmin
  }
  
  else {
    if(max(y) < 0) {
      xmin <- x[y == max(y)]
      if(length(xmin) > 0)
        xmin <- xmin[1]
      xzero <- xmin
    }
    else {
      y1 <- min(y[y > 0])
      if(length(y1) > 0)
        y1 <- y1[length(y1)]
      y2 <- max(y[y < 0])
      if(length(y2) > 0)
        y2 <- y2[1]
      x1 <- x[y == y1]
      if(length(x1) > 0)
        x1 <- x1[length(x1)]
      x2 <- x[y == y2]
      if(length(x2) > 0)
        x2 <- x2[1]
      xzero <- (x2 * y1 - x1 * y2)/(y1 - y2)
      xmin <- x1
      if(abs(y2) < y1)
        xmin <- x2
    }
  }
  resu <-  xzero
  resu
}

#QRLM function for finding area-specific tuning parameters----------------------
QRLM <-function(x,
               y,
               case.weights = rep(1, nrow(x)),
               k=1.345,
               var.weights = rep(1, nrow(x)),
               ..., 
               w = rep(1, nrow(x)),
               init = "ls",
               psi = psi.huber,
               scale.est = c("MAD", "Huber", "proposal 2"),
               k2 = 1.345,
               method = c("M", "MM"),
               maxit = 20,
               acc = 1e-04,
               test.vec = "resid",
               q = 0.5){
  require(MASS)
  irls.delta <- function(old, new) sqrt(sum((old - new)^2)/max(1e-20, sum(old^2)))
  irls.rrxwr <- function(x, w, r) {
    w <- sqrt(w)
    max(abs((matrix(r*w, 1, length(r)) %*% x)/sqrt(matrix(w, 1, length(r)) %*% (x^2))))/sqrt(sum(w*r^2))
  }
  method <- match.arg(method)
  nmx <- deparse(substitute(x))
  if (is.null(dim(x))) {
    x <- as.matrix(x)
    colnames(x) <- nmx
  }
  else x <- as.matrix(x)
  if (is.null(colnames(x)))
    colnames(x) <- paste("X", seq(ncol(x)), sep = "")
  if (qr(x)$rank < ncol(x))
    stop("x is singular: singular fits are not implemented in rlm")
  if (!(any(test.vec == c("resid", "coef", "w", "NULL")) || is.null(test.vec)))
    stop("invalid testvec")
  if (length(var.weights) != nrow(x))
    stop("Length of var.weights must equal number of observations")
  if (any(var.weights < 0))
    stop("Negative var.weights value")
  if (length(case.weights) != nrow(x))
    stop("Length of case.weights must equal number of observations")
  w <- (w * case.weights)/var.weights
  if (method == "M") {
    scale.est <- match.arg(scale.est)
    if (!is.function(psi))
      psi <- get(psi, mode = "function")
    arguments <- list(...)
    if (length(arguments)) {
      pm <- pmatch(names(arguments), names(formals(psi)), nomatch = 0)
      if (any(pm == 0))
        warning(paste("some of ... do not match"))
      pm <- names(arguments)[pm > 0]
      formals(psi)[pm] <- unlist(arguments[pm])
    }
    if (is.character(init)) {
      if (init == "ls")
        temp <- lm.wfit(x, y, w, method = "qr")
      else if (init == "lts")
        temp <- lqs.default(x, y, intercept = FALSE, nsamp = 200)
      else stop("init method is unknown")
      coef <- temp$coef
      resid <- temp$resid
    }
    else {
      if (is.list(init))
        coef <- init$coef
      else coef <- init
      resid <- y - x %*% coef
    }
  }
  else if (method == "MM") {
    scale.est <- "MM"
    temp <- lqs.default(x, y, intercept = FALSE, method = "S", k0 = 1.548)
    coef <- temp$coef
    resid <- temp$resid
    psi <- psi.bisquare
    if (length(arguments <- list(...)))
      if (match("c", names(arguments), nomatch = FALSE)) {
        c0 <- arguments$c
        if (c0 > 1.548) {
          psi$c <- c0
        }
        else warning("c must be at least 1.548 and has been ignored")
      }
    scale <- temp$scale
  }
  else stop("method is unknown")
  done <- FALSE
  conv <- NULL
  n1 <- nrow(x) - ncol(x)
  if (scale.est != "MM")
    scale <- mad(resid/sqrt(var.weights), 0)
  qest <- matrix(0, nrow = ncol(x), ncol = length(q))
  qwt <- matrix(0, nrow = nrow(x), ncol = length(q))
  qfit <- matrix(0, nrow = nrow(x), ncol = length(q))
  qres <- matrix(0, nrow = nrow(x), ncol = length(q))
  qvar.matrix <- array(rep(0,ncol(x)*ncol(x)),dim=c(ncol(x),ncol(x),length(q)))
  qscale <- NULL
  for(i in 1:length(q)) {
    for (iiter in 1:maxit) {
      if (!is.null(test.vec))
        testpv <- get(test.vec)
      if (scale.est != "MM") {
        if (scale.est == "MAD")
          scale <- median(abs(resid/sqrt(var.weights)))/0.6745
        else {gamma<- 4*k2^2*(1-pnorm(k2))*((1-q[i])^2+q[i]^2) - 4*k2*dnorm(k2)*((1-q[i])^2+q[i]^2) + 4*(1-q[i])^2*(pnorm(0)-(1-pnorm(k2))) + 4*q[i]^2*(pnorm(k2)-pnorm(0))
        scale <- sqrt(sum(pmin(resid^2/var.weights,(k2*scale)^2))/(n1*gamma))
        }
        if (scale == 0) {
          done <- TRUE
          break
        }
      }
      w <- psi(resid/(scale * sqrt(var.weights)),k=k) * case.weights
      ww <- 2 * (1 - q[i]) * w
      ww[resid > 0] <- 2 * q[i] * w[resid > 0]
      w <- ww
      temp <- lm.wfit(x, y, w, method = "qr")
      coef <- temp$coef
      resid <- temp$residuals
      if (!is.null(test.vec))
        convi <- irls.delta(testpv, get(test.vec))
      else convi <- irls.rrxwr(x, wmod, resid)
      conv <- c(conv, convi)
      done <- (convi <= acc)
      if (done)
        break
    }
    if (!done)
      warning(paste("rlm failed to converge in", maxit, "steps at q = ", q[i]))
    qest[, i] <- coef
    qscale[i] <-scale
    qwt[, i] <- w
    qfit[, i] <- temp$fitted.values
    qres[,i] <- resid
    
    tmp.res.mq <- qres[,i]/qscale[i]
    Epsi2 <-(sum((qwt[,i]*tmp.res.mq)^2)/(nrow(x)-ncol(x)))
    Epsi <-(1/qscale[i])*(sum(2*(q[i]*(0<=tmp.res.mq & tmp.res.mq<= k)+(1-q[i])*(-k <=tmp.res.mq & tmp.res.mq<0)))/nrow(x))
    qvar.matrix[, , i] <- (((Epsi2)/Epsi^2)*solve(t(x)%*%x))
  }
  list(fitted.values = qfit,
       residuals = qres,
       q.values = q,
       q.weights = qwt,
       coef = qest,
       qscale = qscale,
       var.beta = qvar.matrix)
}

# Original algorithm on Lahiri & Salvati (2023)
# mqre_g<-function(qtl,y,x,group,tol=1e-06,maxit=100,k=100,k_sigma_u=100,k_sigma_e=100)
# {
#   ###############functions################
#   
#   my.psi<-function(u,q,k){
#     sm<-median(abs(u))/0.6745
#     w <- psi.huber(u/sm,k)
#     ww <- 2 * (1 - q) * w
#     ww[u> 0] <- 2 * q * w[u > 0]
#     w <- ww
#     w*u
#   }
#   
#   der.psi<-function(u,q,k){
#     sm<-median(abs(u))/0.6745
#     u<-u/sm
#     der.psi <- ifelse(abs(u) <= k, 1, 0)
#     ww <- 2 * (1 - q) * der.psi
#     ww[u> 0] <- 2 * q * der.psi[u > 0]
#     w <- ww
#     w
#   }
#   ########################starting data###################
#   n=length(y)
#   ni=table(group)
#   m=length(ni)
#   areanumber=m
#   p=ncol(x)
#   gg=unique(group)
#   z<-model.matrix(y~as.factor(group)-1)
#   
#   #STEP 1 (Definition of the starting values)
#   w=cbind(x,z)
#   fit.H=lm(y~w)
#   e.1.H=residuals(fit.H)
#   sigma.e.hat.H=sum(e.1.H^2)/(n-(qr(w)$rank))
#   
#   xx=x
#   fit2.H=lm(y~xx)
#   e.2.H=residuals(fit2.H)
#   A.H= sum(e.2.H^2)
#   xx=as.matrix(xx)
#   B.H=sigma.e.hat.H*(n-qr(xx)$rank)
#   o.H=diag(n)-(xx%*%(ginv(t(xx)%*%xx)%*%t(xx)))
#   C.H=sum(diag(t(z)%*%o.H%*%z))
#   sigma.v.hat.H=(A.H-B.H)/C.H
#   sigmasq0= sigma.e.hat.H #initial values of sigma_sq_e#
#   sigmasq0.v=sigma.v.hat.H #initial values sigma_sq_V#
#   
#   if(dim(x)[2]==1){ls <- lm(y ~ x[,1]-1)}else{
#     ls <- lm(y ~ x[,-1])}
#   beta1<-c(ls$coefficients)
#   estsigma2u<-sigmasq0.v
#   estsigma2e<-sigmasq0
#   sigma1G.old<-NULL
#   sigma1E.old<-NULL
#   sigma1E.old<-rep(sigmasq0,areanumber)
#   sigma1G.old<-c(sigmasq0.v)
#   sigma1G.est<-NULL
#   sigma1E.est<-NULL
#   
#   iter<-0
#   ZZ<-z%*%t(z)
#   beta.old<-matrix(beta1,areanumber,p,byrow=T)
#   beta.est<-matrix(0,areanumber,p)
#   var.beta<-array(rep(0,ncol(x)*ncol(x)),dim=c(ncol(x),ncol(x),length(qtl)))
#   
#   #STEP 2: estimation of betas and sigma2
#   while(TRUE)
#   { residuals<-NULL
#   for (i in 1:length(qtl))
#   {
#     #STEP 1 Estimation of beta
#     xbeta <- c(x %*% beta.old[i,])
#     R <- diag(rep(sigma1E.old[i], n),n,n)
#     G <- diag(rep(sigma1G.old, areanumber),areanumber,areanumber)
#     V <- R + z %*% G %*% t(z)
#     V.inv<-chol2inv(chol(V))
#     U <- diag(diag(V))
#     U.inv <- chol2inv(chol(U))
#     ystar<-c(as.numeric(sqrt(U.inv) %*% y))
#     xstar<-matrix(as.numeric(sqrt(U.inv) %*% x),n,p)
#     mod<-QRLM(y=ystar,x=xstar,q=qtl[i],maxit=maxit,acc=tol,k=k)
#     beta.est[i,]<-mod$coef
#     var.beta[,,i]<-mod$var.beta
#     r <- c(sqrt(U.inv) %*% (y - xbeta))
#     #mod.rlmerE<-rlmer(r~1+(1|group),rho.sigma.e = psi2propII(smoothPsi, k = k_sigma),
#     #   rho.sigma.b = psi2propII(smoothPsi, k = k_sigma))
#     psi.r <- my.psi(as.vector(r), k = k_sigma_e, q=qtl[i])
#     qq <- V.inv %*% sqrt(U) %*% psi.r
#     KK5<-sqrt(U) %*% V.inv
#     KK0<-t(psi.r)%*% KK5
#     A1=KK0%*%qq
#     A2=KK0%*%ZZ%*%qq
#     A <- matrix(c(A1[1],A2[1]), nrow = 2, ncol=1)
#     const<- 4*k_sigma_e^2*(1-pnorm(k_sigma_e))*((1-qtl[i])^2+qtl[i]^2) - 4*k_sigma_e*dnorm(k_sigma_e)*((1-qtl[i])^2+qtl[i]^2) + 4*(1-qtl[i])^2*(pnorm(0)-(1-pnorm(k_sigma_e))) + 4*qtl[i]^2*(pnorm(k_sigma_e)-pnorm(0))
#     K2_value<- const
#     K2<-diag(K2_value, n)
#     KK1<-K2%*%V.inv
#     KK2<-KK1%*%V.inv
#     KK3<-KK1%*%ZZ%*%V.inv
#     t1=sum(diag(KK2))
#     t2=sum(diag(KK2%*%ZZ))
#     t3=sum(diag(KK3))
#     t4=sum(diag(KK3%*%ZZ))
#     T1<- matrix(c(t1,t3,t2,t4), nrow = 2, ncol=2)
#     sigma=chol2inv(chol(T1))%*%A
#     sigma1E.est[i]<-sigma[1,1]
#     if(p==1)residuals<-c(residuals,(y[group==i]-x[group==i,1]*beta.est[i,1]))
#     if(p==2)residuals<-c(residuals,(y[group==i]-x[group==i,2]*beta.est[i,2]))
#     if(p>2)residuals<-c(residuals,(y[group==i]-x[group==i,-1]%*%beta.est[i,-1]))
#   }#end beta estimation
#   
#   #STEP 3: computation variance components sigma^2_gamma
#   residuals<-residuals-mean(residuals)
#   R <- diag(rep(sigma1E.old, ni),n,n)
#   G <- diag(rep(sigma1G.old, areanumber),areanumber,areanumber)
#   V <- R + z %*% G %*% t(z)
#   V.inv<-chol2inv(chol(V))
#   U <- diag(diag(V))
#   U.inv <- chol2inv(chol(U))
#   r <- c(sqrt(U.inv) %*% residuals)
#   psi.r <- my.psi(as.vector(r), k = k_sigma_u, q=0.5)
#   qq <- V.inv %*% sqrt(U) %*% psi.r
#   KK5<-sqrt(U) %*% V.inv
#   KK0<-t(psi.r)%*% KK5
#   A1=KK0%*%qq
#   A2=KK0%*%ZZ%*%qq
#   A <- matrix(c(A1[1],A2[1]), nrow = 2, ncol=1)
#   const<- 4*k_sigma_u^2*(1-pnorm(k_sigma_u))*((1-0.5)^2+0.5^2) - 4*k_sigma_u*dnorm(k_sigma_u)*((1-0.5)^2+0.5^2) + 4*(1-0.5)^2*(pnorm(0)-(1-pnorm(k_sigma_u))) + 4*0.5^2*(pnorm(k_sigma_u)-pnorm(0))
#   K2_value<- const
#   K2<-diag(K2_value, n)
#   KK1<-K2%*%V.inv
#   KK2<-KK1%*%V.inv
#   KK3<-KK1%*%ZZ%*%V.inv
#   t1=sum(diag(KK2))
#   t2=sum(diag(KK2%*%ZZ))
#   t3=sum(diag(KK3))
#   t4=sum(diag(KK3%*%ZZ))
#   T1<- matrix(c(t1,t3,t2,t4), nrow = 2, ncol=2)
#   sigma=chol2inv(chol(T1))%*%A
#   sigma1G.est<-sigma[2,1]
#   if(p==1)B0<-0
#   if(p>1)B0<-mean(beta.est[,1])
#   
#   if(p>1)fehler <- sum((beta.old[,-1] -beta.est[,-1])^2)+ sum((sigma1G.old -sigma1G.est)^2+sum((sigma1E.old -sigma1E.est)^2))
#   if(p==1)fehler <- sum((beta.old[,1] -beta.est[,1])^2)+ sum((sigma1G.old -sigma1G.est)^2+sum((sigma1E.old -sigma1E.est)^2))
#   ifelse (iter < maxit, ifelse(fehler > tol, {iter=iter+1}, {break}), {break})
#   beta.old <- beta.est
#   sigma1G.old<-sigma1G.est
#   sigma1E.old<-sigma1E.est
#   } #end while, estimation procedure
#   
#   if(sigma1G.est<0){lme_fit<-lmer(y~-1+x+(1|group))
#   sigma1G.est<-summary(lme_fit)$varcor$group[1]}
#   if(p>1)coefficients =cbind(rep(B0,areanumber),beta.est[,-1])
#   if(p==1)coefficients =beta.est
#   list(coefficients =coefficients,sigma2u=(sigma1G.est),sigma2e=(sigma1E.est),varBeta=var.beta,quantile=qtl)
# }

#Function for obtaining beta_i, sigma2e_i and sigma2u---------------------------
mqre_g <- function(fixed, framework, transformation_par, Qi, var.weights, maxit, acc, k_b, k_sigma_e, k_sigma_u){
  model.f <- model.frame(fixed, data = transformation_par$transformed_data)
  ys <- as.numeric(model.response(model.f))
  xs <- model.matrix(fixed, data=transformation_par$transformed_data)[, , drop = FALSE]
  area.s <- framework$smp_domains_vec
  area_id <- unique(area.s)
  #psi-q function
  psi.q<-function(u, q, k){
    sm<-median(abs(u))/0.6745
    w <- psi.huber(u/sm, k)
    ww <- 2 * (1 - q) * w
    ww[u> 0] <- 2 * q * w[u > 0]
    w <- ww
    w*u
  }
  # Obtain beta_i for in-sample domains
  out1 <- QRLM(y = ys, x = xs, q = Qi, var.weights = var.weights, maxit = maxit, acc = acc, k = k_b)

  # Obtain sigmae2_i for in-sample domains
  est.sigma2e <- rep(1, framework$N_dom_smp)
  for (i in 1:framework$N_dom_smp)
  {
    e.1.H <- psi.q(out1$residuals[,i], Qi[i], k = k_sigma_e)
    mod.lmm <- summary(lmer(e.1.H ~ 1 + (1|area.s)))
    est.sigma2e[i] <- as.numeric(mod.lmm$sigma^2)
  }
  # Obtain common sigmau2 for all the domains
  resid <- NULL
  if(ncol(xs) == 2){
    for (i in 1:framework$N_dom_smp)
    {
      pred <- xs[area.s == area_id[i], -1] * out1$coef[-1, i]
      resid <- c(resid, c(ys[area.s == area_id[i]] - pred))
    }
  }
  if(ncol(xs) > 2){
    for (i in 1:framework$N_dom_smp)
    {
      pred <- xs[area.s == area_id[i], -1] %*% out1$coef[-1, i]
      resid <- c(resid, c(ys[area.s == area_id[i]] - pred))
    }
  }
  
  require(emdi)
  res.fh <- tapply(resid, area.s, mean)
  dta.fh <- data.frame(res.fh = as.numeric(res.fh), vardir.fh = as.numeric(est.sigma2e/framework$n_smp), area.fh = area_id)
  fh_std <- fh(fixed = res.fh ~ 1, vardir = "vardir.fh", domains = "area.fh", method = "reblup", k=k_sigma_u, combined_data = dta.fh)
  B0.coef <- as.numeric(fh_std$model$coefficients[1])
  est.sigma2u <- as.numeric(fh_std$model$variance)
  
  # mod.lmm1 <- summary(lmer(resid ~ 1 + (1|area.s)))
  # B0.coef <- mod.lmm1$coefficients[1]
  # est.sigma2u <- as.numeric(mod.lmm1$varcor$area.s[1])
  coefficients <- rbind(rep(B0.coef, framework$N_dom_smp), out1$coef[-1, ])
  colnames(coefficients) <- c(unique(framework$smp_domains_vec))
  
  return(list(coefficients = coefficients,
              sigma2u = (est.sigma2u),
              sigma2e = (est.sigma2e),
              quantile = Qi))
}
#Function for extracting the estimated parameters-------------------------------
model_par_hdp <- function(framework,
                      transformation_par,
                      mixed_model,
                      fixed, 
                      betas.out){
  #Extract response variable and covariates
  model.f <- model.frame(fixed, data = transformation_par$transformed_data)
  ys <- as.numeric(model.response(model.f))
  xs <- model.matrix(fixed, data=transformation_par$transformed_data)[, , drop = FALSE]
  # fixed parameters for in-sample domains
  betas <- mixed_model$coefficients
  # fixed parameters for out-sample domains
  betas.out <- c(betas[1,1], betas.out[-1,])
  # Estimated sampling variance
  sigmae2est <- mixed_model$sigma2e
  # Estimated variance of random effect
  sigmau2est <- mixed_model$sigma2u
  # Random effect: vector with zeros for all domains, filled with 0
  rand_eff <- rep(0, length(unique(framework$pop_domains_vec)))
  # In-sample-domain means for the response variable
  bary <- as.vector(unlist(lapply(split(ys, framework$smp_domains_vec), mean)))
  # In-sample-domain means for covariates
  barx <- matrix(unlist(lapply(split.data.frame(xs, framework$smp_domains_vec), colMeans)),
               framework$N_dom_smp, ncol(xs), byrow = T)
  Bi <- (sigmae2est / framework$n_smp) / (sigmae2est / framework$n_smp + c(sigmau2est))
  rand_eff[framework$dist_obs_dom] <- (1-Bi) * (bary - diag(barx%*%betas))
  
  return(list(
    betas = betas,
    betas.out = betas.out,
    sigmae2est = sigmae2est,
    sigmau2est = sigmau2est,
    rand_eff = rand_eff
  ))
}

# Function gen_model calculates the parameters in the generating model.---------
gen_model_hdp <- function(fixed,
                      framework,
                      model_par) {
  # Parameter for calculating variance of new random effect
  gamma <- model_par$sigmau2est / (model_par$sigmau2est +
                                     model_par$sigmae2est / framework$n_smp)
  # Variance of new random effect
  sigmav2est <- model_par$sigmau2est * (1 - gamma)
  # Random effect in constant part of y for in-sample households
  rand_eff_pop <- rep(model_par$rand_eff, framework$n_pop)
  # Model matrix for population covariate information
  framework$pop_data[[paste0(fixed[2])]] <- seq_len(nrow(framework$pop_data))
  X_pop <- model.matrix(fixed, framework$pop_data)
  
  # Constant part of predicted y
  mu_fixed <- vector(length = framework$N_pop)
  # synthetic part for in-sample-domains
  betas_in <- apply(model_par$betas, 1, rep, framework$n_pop[framework$dist_obs_dom])
  mu_fixed[framework$obs_dom] <- diag(X_pop[framework$obs_dom,] %*% t(betas_in))
  # synthetic part for out-sample-domains
  mu_fixed[!framework$obs_dom] <- X_pop[!framework$obs_dom,] %*% model_par$betas.out
  
  mu <- mu_fixed + rand_eff_pop
  
  return(list(sigmav2est = sigmav2est,
              mu = mu,
              mu_fixed = mu_fixed))
}

# The function errors_gen_hdp returns error terms of the generating model.----------
errors_gen_hdp <- function(framework, model_par, gen_model){
  epsilon <- rep(0, framework$N_pop)
  epsilon[framework$obs_dom] <- rnorm(sum(framework$n_pop[framework$dist_obs_dom]),
                                      0, 
                                      rep(sqrt(model_par$sigmae2est), framework$n_pop[framework$dist_obs_dom]))
  vu <- vector(length = framework$N_pop)
  # new random effect for out-of-sample domains
  vu[!framework$obs_dom] <- rep(
    rnorm(
      framework$N_dom_unobs,
      0,
      sqrt(model_par$sigmau2est)
    ),
    framework$n_pop[!framework$dist_obs_dom]
  )
  # new random effect for in-sample-domains
  vu[framework$obs_dom] <- rep(
    rnorm(
      rep(1, framework$N_dom_smp),
      0,
      sqrt(gen_model$sigmav2est)
    ),
    framework$n_pop[framework$dist_obs_dom]
  )
  # individual error term in generating model epsilon
  return(list(epsilon = epsilon,
              vu = vu))
}

# The function prediction_y returns a predicted income vector which can be used
# to calculate indicators. Note that a whole income vector is predicted without
# distinction between in- and out-of-sample domains.
prediction_y <- function(gen_model,
                         errors_gen,
                         framework,
                         fixed,
                         transformation,
                         lambda,
                         shift) {
  
  # predicted population income vector
  y_pred <- gen_model$mu + errors_gen$epsilon + errors_gen$vu
  
  # back-transformation of predicted population income vector
  y_pred <- back_transformation(
    y = y_pred,
    transformation = transformation,
    lambda = lambda,
    shift = shift,
    framework = framework,
    fixed = fixed
  )
  y_pred[!is.finite(y_pred)] <- 0
  
  return(y_pred)
} # End prediction_y

# Monte-Carlo approximation ----------------------------------------------------
monte_carlo <- function(transformation,
                        L,
                        framework,
                        lambda,
                        shift,
                        model_par,
                        gen_model,
                        fixed) {
  N_dom_pop_tmp <- framework$N_dom_pop
  pop_domains_vec_tmp <- framework$pop_domains_vec
  
  ests_mcmc <- array(dim = c(
    N_dom_pop_tmp,
    L,
    length(framework$indicator_names)
  ))
  
  for (l in seq_len(L)) {
    # Errors in generating model: individual error term and random effect
    # See below for function errors_gen.
    errors <- errors_gen_hdp(
      framework = framework,
      model_par = model_par,
      gen_model = gen_model
    )
    
    # Prediction of population vector y
    # See below for function prediction_y.
    population_vector <- prediction_y(
      gen_model = gen_model,
      errors_gen = errors,
      framework = framework,
      fixed = fixed,
      transformation = transformation,
      lambda = lambda,
      shift = shift
    )
    
    if(!is.null(framework$pop_weights)){
      pop_weights_vec <- framework$pop_data[[framework$pop_weights]]
    }else{
      pop_weights_vec <- rep(1, nrow(framework$pop_data))
    }
    
    # Calculation of indicators for each Monte Carlo population
    ests_mcmc[, l, ] <-
      matrix(
        nrow = N_dom_pop_tmp,
        data = unlist(lapply(framework$indicator_list,
                             function(f, threshold) {
                               matrix(
                                 nrow = N_dom_pop_tmp,
                                 data = unlist(mapply(
                                   y = split(population_vector, pop_domains_vec_tmp),
                                   pop_weights = split(pop_weights_vec, pop_domains_vec_tmp),
                                   f,
                                   threshold = framework$threshold
                                 )), byrow = TRUE
                               )
                             },
                             threshold = framework$threshold
        ))
      )
  } # End for loop
    
  # Point estimations of indicators by taking the mean
  point_estimates <- data.frame(
    Domain = unique(pop_domains_vec_tmp),
    apply(ests_mcmc, c(3), rowMeans)
  )
  colnames(point_estimates) <- c("Domain", framework$indicator_names)
  return(point_estimates)
} # End Monte-Carlo

