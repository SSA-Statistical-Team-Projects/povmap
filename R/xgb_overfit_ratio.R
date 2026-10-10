#' Out-of-fold to in-sample residual ratio of tuning configurations
#'
#' For each candidate configuration, compares how well the fitted booster
#' reproduces the sample it was trained on with how well it predicts clusters
#' it has not seen. A configuration that fits the sample far more closely than
#' it predicts held-out clusters is overfit; the ratio is the measure used by
#' \code{xgb_top_configs(rule = "least_overfit")}.
#'
#' For every configuration (one row of \code{configs}) the booster is fitted
#' once on the whole sample and the in-sample residuals are
#' \code{y - fitted}, on the transformed scale. The out-of-fold residuals come
#' from \code{folds}-fold cross-validation with whole clusters held out, the
#' fold of each cluster drawn exactly as \code{\link{xgb_tune}} draws it
#' (\code{set.seed(seed)}, then one fold for each unique cluster in order of
#' appearance). Every model, the full one and the fold ones, is weighted as
#' \code{\link{xgb}} would weight a fit on its own rows (see
#' \code{rescale_weights} and \code{variance_y} in \code{\link{xgb_tune}}).
#' The residual spread is the survey-weighted standard deviation, in three
#' pools: \code{total} (all observations), \code{domain} (the residuals
#' averaged within each domain, weighted by the domain's weight total) and
#' \code{sub} (the residuals less their domain mean, the within-domain part).
#' The ratio of each is the out-of-fold spread divided by the in-sample one.
#'
#' This is the procedure of the DRC \dQuote{least overfit} selection. There the
#' ratios were computed in a separate pass on a fresh fold assignment (seed 901)
#' rather than on the folds that chose the candidates, with xgboost fit seed 123
#' and one thread per fit, which are the defaults here for \code{fit_seed} and
#' the threads. The scores used for the one-SE set came from two such fresh
#' assignments (seeds 901 and 902); see \code{\link{xgb_top_configs}}.
#'
#' @inheritParams xgb_tune
#' @param configs data.frame of configurations, one per row, with the
#'   hyperparameter columns of \code{\link{xgb_tune}} (\code{nround} or
#'   \code{nrounds}, \code{max_depth}, \code{colsample_bytree},
#'   \code{colsample_bylevel}, \code{colsample_bynode}, \code{subsample},
#'   \code{min_child_weight}, \code{eta}, \code{gamma}, \code{max_delta_step},
#'   \code{lambda}, \code{alpha}), for example the \code{candidates} of an
#'   \code{xgb_tune} result. Missing hyperparameters take \code{xgb()}'s
#'   defaults.
#' @param transformation \code{"no"}, \code{"arcsin"}, \code{"log"} or
#'   \code{"poisson"} (log1p), the scale on which the booster is trained and
#'   the residuals are taken. Defaults to \code{"no"}.
#' @param rows optional integer vector of the rows of \code{configs} to
#'   evaluate (default all). Rows not evaluated are not returned.
#' @param seed integer seed for the assignment of clusters to folds, as in
#'   \code{\link{xgb_tune}}. \code{NULL} leaves the random number stream as it is.
#' @param fit_seed seed passed to xgboost for every fit, so that the
#'   configurations share common random numbers. Default 123. \code{NULL}
#'   leaves xgboost's own default.
#' @param cpus number of worker processes (PSOCK), default 1.
#' @param ... further parameters passed to \code{xgb.train}.
#' @return a data.frame with one row per evaluated configuration:
#'   \code{row} (its row in \code{configs}), \code{in_total}, \code{oof_total},
#'   \code{in_domain}, \code{oof_domain}, \code{in_sub}, \code{oof_sub} (the
#'   in-sample and out-of-fold residual standard deviations) and
#'   \code{ratio_total}, \code{ratio_domain}, \code{ratio_sub}
#'   (out-of-fold over in-sample).
#' @seealso \code{\link{xgb_top_configs}}, \code{\link{xgb_tune}}
#' @export
#' @importFrom foreach foreach %dopar%
#' @importFrom doParallel registerDoParallel
xgb_overfit_ratio <- function(fixed, smp_data, configs, smp_weights = NULL, domains,
                              cluster = "domains", transformation = "no", folds = 10,
                              rescale_weights = TRUE, variance_y = NULL, rows = NULL,
                              seed = NULL, fit_seed = 123, cpus = 1, ...) {
  .assert_xgb_version()
  outcome <- all.vars(fixed[[2]])
  covariates <- all.vars(fixed[[3]])
  if (identical(cluster, "domains")) cluster <- domains
  X_smp <- smp_data[, unique(c(covariates, domains, cluster)), drop = FALSE]
  Y_smp <- data.frame(labels = smp_data[[outcome]])
  w <- if (is.null(smp_weights)) rep(1, nrow(Y_smp)) else as.numeric(smp_data[[smp_weights]])
  xgb_check2(transformation = transformation, Y_smp = Y_smp, X_smp = X_smp,
             smp_weights = w, domains = domains, cluster = cluster)
  yt <- .xgb_transform_y(Y_smp$labels, transformation)
  het <- NULL
  if (!is.null(variance_y)) {
    if (!is.character(variance_y) || length(variance_y) != 1L || !variance_y %in% names(smp_data))
      stop("variance_y must name one variable in smp_data.", call. = FALSE)
    het <- as.numeric(smp_data[[variance_y]])^-0.5
    if (anyNA(het) || any(!is.finite(het)))
      stop("variance_y must be numeric, positive and not missing.", call. = FALSE)
  }
  grid <- .xgb_ratio_grid(configs)
  cl_ids <- X_smp[[cluster]]
  if (!is.null(seed)) set.seed(seed)
  cu <- unique(cl_ids)
  fold_of <- sample(seq_len(folds), length(cu), replace = TRUE)
  fold <- fold_of[match(cl_ids, cu)]
  .xgb_overfit_core(grid, X = X_smp[, setdiff(covariates, c(domains, cluster)), drop = FALSE],
                    yt = yt, w = w, dom = X_smp[[domains]], fold = fold, folds = folds,
                    rescale_weights = rescale_weights, het = het, fit_seed = fit_seed,
                    rows = rows, cpus = cpus, dots = list(...))
}

## transformed outcome, as xgb_tune() and xgb() train on it
.xgb_transform_y <- function(y, transformation) {
  switch(transformation,
         no = y, arcsin = asin(sqrt(y)), log = log(y), poisson = log1p(y),
         stop("transformation must be one of no, arcsin, log or poisson for the overfit ratio.", call. = FALSE))
}

## hyperparameter table with xgb_tune's names, missing ones from xgb()'s defaults
.xgb_ratio_grid <- function(configs) {
  if (!is.data.frame(configs) || nrow(configs) < 1)
    stop("configs must be a data.frame with one configuration per row.", call. = FALSE)
  if ("nrounds" %in% names(configs) && !"nround" %in% names(configs))
    names(configs)[names(configs) == "nrounds"] <- "nround"
  d <- .xgb_xgb_defaults(); names(d)[names(d) == "nrounds"] <- "nround"
  out <- data.frame(row.names = seq_len(nrow(configs)))
  for (h in .xgb_hp_names) {
    v <- if (h %in% names(configs)) configs[[h]] else rep(d[[h]], nrow(configs))
    if (!is.numeric(v) || anyNA(v)) stop("configs$", h, " must be numeric without NA.", call. = FALSE)
    out[[h]] <- as.numeric(v)
  }
  out
}

## weighted standard deviation (population form, as the DRC selection computes it)
.xgb_wsd <- function(x, w) { mu <- sum(w * x) / sum(w); sqrt(sum(w * (x - mu)^2) / sum(w)) }

## Spread of residuals `res` in the total, domain and within-domain pools
.xgb_resid_pools <- function(res, w, dom) {
  dm <- stats::ave(res * w, dom, FUN = sum) / stats::ave(w, dom, FUN = sum)
  tw <- tapply(w, dom, sum)
  tm <- tapply(res * w, dom, sum) / tw
  c(total = .xgb_wsd(res, w), domain = .xgb_wsd(as.numeric(tm), as.numeric(tw)),
    sub = .xgb_wsd(res - dm, w))
}

## One configuration: the full fit and the cluster-grouped fold fits
.xgb_overfit_one <- function(g, X, yt, w, dom, fold, folds, rescale_weights, het, fit_seed, dots) {
  p <- c(list(max_depth = g$max_depth, colsample_bytree = g$colsample_bytree,
              colsample_bylevel = g$colsample_bylevel, colsample_bynode = g$colsample_bynode,
              subsample = g$subsample, min_child_weight = g$min_child_weight, eta = g$eta,
              gamma = g$gamma, max_delta_step = g$max_delta_step, lambda = g$lambda,
              alpha = g$alpha, nthread = 1), dots)
  if (!is.null(fit_seed)) p$seed <- fit_seed
  XM <- data.matrix(X)
  fitm <- function(tr) xgboost::xgb.train(
    data = xgboost::xgb.DMatrix(XM[tr, , drop = FALSE], label = yt[tr],
                                weight = .xgb_fit_weights(w[tr], dom[tr], rescale_weights, het[tr])),
    params = p, nrounds = g$nround, verbose = 0)
  r_in <- yt - predict(fitm(rep(TRUE, length(yt))), XM)
  pr <- numeric(length(yt))
  for (k in seq_len(folds)) {
    h <- fold == k
    if (any(h)) pr[h] <- predict(fitm(!h), XM[h, , drop = FALSE])
  }
  a <- .xgb_resid_pools(r_in, w, dom); b <- .xgb_resid_pools(yt - pr, w, dom)
  c(in_total = a[["total"]], oof_total = b[["total"]], in_domain = a[["domain"]],
    oof_domain = b[["domain"]], in_sub = a[["sub"]], oof_sub = b[["sub"]])
}

## Ratios for the rows `rows` of `grid` (xgb_tune names); `fold` gives each observation's fold
.xgb_overfit_core <- function(grid, X, yt, w, dom, fold, folds, rescale_weights = TRUE, het = NULL,
                              fit_seed = 123, rows = NULL, cpus = 1, dots = list()) {
  if (is.null(rows)) rows <- seq_len(nrow(grid))
  if (is.null(het)) het <- rep(1, length(yt))
  one <- .xgb_overfit_one
  cl <- parallel::makeCluster(cpus)
  doParallel::registerDoParallel(cl)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  res <- foreach::foreach(r = rows, .combine = rbind, .packages = "xgboost") %dopar% {
    c(row = r, one(grid[r, , drop = FALSE], X, yt, w, dom, fold, folds, rescale_weights, het,
                   fit_seed, dots))
  }
  res <- matrix(as.numeric(res), ncol = 7)
  d <- as.data.frame(res)
  names(d) <- c("row", "in_total", "oof_total", "in_domain", "oof_domain", "in_sub", "oof_sub")
  d$row <- as.integer(d$row)
  d$ratio_total <- d$oof_total / d$in_total
  d$ratio_domain <- d$oof_domain / d$in_domain
  d$ratio_sub <- d$oof_sub / d$in_sub
  d
}
