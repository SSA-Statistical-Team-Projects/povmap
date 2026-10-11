# Out-of-fold residual pools for xgb() (area_pool, sub_pool, *_deflate, pool_tune).
#
# xgb() generates the bootstrap population from two residual pools: an area pool (one value per
# sampled domain: the weighted mean residual of the domain's units) and a sub-area pool (a unit's
# residual less its domain's mean). By default both come from the final fit's own residuals
# ("in_sample"), which are smaller than the errors the model makes on units it has not seen. The
# "out_of_fold" option takes them from models that did not train on the unit. Everything here
# runs AFTER the point estimate, so point estimates do not depend on the pool choice, and it saves
# and restores the random-number state, so the bootstrap streams are those of the default call.

## Fold of every sample row, domain-grouped. `assignment` (optional): a data.frame whose first two
## columns are the domain label and its fold, e.g. the fold_assignment of an xgb_tune result.
.xgb_pool_fold_domain <- function(dom, K, seed, assignment = NULL) {
  dom <- as.character(dom)
  if (!is.null(assignment)) {
    a <- as.data.frame(assignment)
    if (ncol(a) < 2L) stop("pool_folds as a data.frame needs the domain label and its fold in the first two columns.", call. = FALSE)
    f <- as.integer(a[[2]])[match(dom, as.character(a[[1]]))]
    if (anyNA(f)) stop("pool_folds: no fold for domain(s) ", paste(utils::head(unique(dom[is.na(f)]), 5), collapse = ", "), ".", call. = FALSE)
    return(f)
  }
  ud <- unique(dom)
  if (K < 2L || K > length(ud)) stop("pool_folds must be between 2 and the number of domains (", length(ud), ").", call. = FALSE)
  had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  if (had) old <- get(".Random.seed", envir = globalenv())
  set.seed(seed)
  fu <- sample(rep(seq_len(K), length.out = length(ud)))
  if (had) assign(".Random.seed", old, envir = globalenv()) else rm(".Random.seed", envir = globalenv())
  fu[match(dom, ud)]
}

## Fold of every sample row, grouped by sub-area within domains: a domain's sub-areas get distinct folds (cycling
## when a domain has more sub-areas than folds), so a held-out sub-area's fit contains its domain's other sub-areas
## except those that share its fold (at most ceiling(n_sub / K) - 1).
.xgb_pool_fold_sub <- function(cell, dom, K, seed) {
  cell <- as.character(cell); dom <- as.character(dom)
  had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  if (had) old <- get(".Random.seed", envir = globalenv())
  set.seed(seed)
  fold <- integer(length(cell))
  for (d in unique(dom)) {
    cells <- unique(cell[dom == d])
    f <- sample(rep_len(sample.int(K), length(cells)))
    fold[dom == d] <- f[match(cell[dom == d], cells)]
  }
  if (had) assign(".Random.seed", old, envir = globalenv()) else rm(".Random.seed", envir = globalenv())
  fold
}

## Out-of-fold prediction (model scale) of every sample row at its census cell's covariates `Xcell`.
## Each fold model is fitted on the other folds exactly as xgb() fits: weights rescaled within domain on the
## training rows, divided by their mean, times the heteroscedasticity factor.
.xgb_oof_hat <- function(fold, X, y_t, w, dom, het, Xcell, params, nrounds, rescale_weights, seed, fit_seed) {
  hat <- rep(NA_real_, length(y_t)); X <- as.matrix(X)
  had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  if (had) old <- get(".Random.seed", envir = globalenv())
  for (k in sort(unique(fold))) {
    te <- which(fold == k); tr <- which(fold != k)
    p <- params; if (!is.null(fit_seed)) p$seed <- fit_seed
    set.seed(seed)
    dtrain <- xgboost:::xgb.DMatrix(data = X[tr, , drop = FALSE], label = y_t[tr],
                                    weight = .xgb_fit_weights(w[tr], dom[tr], rescale_weights, het[tr]))
    fit <- xgboost:::xgb.train(data = dtrain, params = p, nrounds = nrounds, verbose = 0)
    te <- te[!is.na(Xcell[te, 1])]
    if (length(te)) hat[te] <- predict(fit, Xcell[te, , drop = FALSE])
  }
  if (had) assign(".Random.seed", old, envir = globalenv()) else rm(".Random.seed", envir = globalenv())
  hat
}

## Pools from unit residuals `r` (weights `w`, domain labels `g`), arithmetic as point_estim_xgb()
.xgb_pools_from_resid <- function(r, w, g) {
  rd <- collapse:::fmean(r, g = g, w = w, use.g.names = TRUE)
  wts <- collapse:::fsum(w, g = g)
  list(resid_sub_domains = r - rd[match(as.character(g), names(rd))], resid_domains = as.numeric(rd),
       wt_sub_domains = w, wt_domains = as.numeric(wts), names = names(rd))
}

## Response-scale variance -> model scale, variance * slope^2, with the package's delta-method slope
.xgb_var_to_model_scale <- function(v, est, transformation, transform_outcome) {
  out <- rep(NA_real_, length(v))
  for (i in which(!is.na(v) & !is.na(est)))
    out[i] <- v[i] * .perturb_slope(est[i], sqrt(v[i]), transformation, transform_outcome)^2
  out
}

## Cluster-robust (Taylor linearisation, sample rows as the clusters, with-replacement form) variance of each
## domain's weighted mean: n/(n - 1) * sum((w_i (y_i - ybar))^2) / (sum w)^2; NA for a domain with one row.
.xgb_domain_mean_var <- function(y, w, g) {
  g <- as.character(g); out <- stats::setNames(rep(NA_real_, length(unique(g))), unique(g))
  for (d in names(out)) {
    i <- g == d; n <- sum(i)
    if (n < 2L) next
    yb <- sum(w[i] * y[i]) / sum(w[i])
    out[d] <- n / (n - 1) * sum((w[i] * (y[i] - yb))^2) / sum(w[i])^2
  }
  out
}

## Deflate pool values `r` (weights w) by sqrt(max(msr - msv, 0) / msr); msr and msv are weighted means over the
## entries that have a sampling variance (v_t, model scale).
.xgb_deflate_pool <- function(r, w, v_t, what) {
  av <- !is.na(v_t)
  if (!any(av)) stop("pool deflation (", what, "): no entry of the pool has a sampling variance.", call. = FALSE)
  if (any(!av)) warning("pool deflation (", what, "): ", sum(!av), " of ", length(av),
                        " pool entries have no sampling variance and are left out of the two means.", call. = FALSE)
  msr <- stats::weighted.mean(r[av]^2, w[av]); msv <- stats::weighted.mean(v_t[av], w[av])
  f <- if (msr > 0) sqrt(max(msr - msv, 0) / msr) else 1
  list(resid = r * f, factor = f, msr = msr, msv = msv, n_missing = sum(!av))
}

## sampling variance input -> vector aligned with the given keys; `kind` "area" (by domain) or "sub" (by row)
.xgb_resolve_area_var <- function(spec, fwk, ok, dom_pop, nm, y, w) {
  if (is.null(spec)) return(.xgb_domain_mean_var(y, w, dom_pop)[nm])
  if (is.character(spec) && length(spec) == 1L) {
    if (!spec %in% names(fwk$smp_data)) stop("area_sampling_var: no column '", spec, "' in smp_data.", call. = FALSE)
    col <- fwk$smp_data[[spec]][ok]
    first <- tapply(col, dom_pop, function(z) z[1]); const <- tapply(col, dom_pop, function(z) all(z == z[1] | is.na(z)))
    if (!all(const)) stop("area_sampling_var: column '", spec, "' is not constant within domain.", call. = FALSE)
    return(as.numeric(first[nm]))
  }
  if (!is.numeric(spec) || is.null(names(spec)))
    stop("area_sampling_var must be NULL, a column name of smp_data, or a named numeric vector (names = domain labels).", call. = FALSE)
  as.numeric(spec[match(nm, names(spec))])
}
.xgb_resolve_sub_var <- function(spec, fwk, ok, variance_y) {
  if (is.null(spec)) {
    if (is.null(variance_y)) stop("sub_deflate needs sub_sampling_var (or variance_y, whose column is then used).", call. = FALSE)
    spec <- variance_y
  }
  v <- if (is.character(spec) && length(spec) == 1L) {
    if (!spec %in% names(fwk$smp_data)) stop("sub_sampling_var: no column '", spec, "' in smp_data.", call. = FALSE)
    fwk$smp_data[[spec]]
  } else spec
  if (!is.numeric(v) || length(v) != nrow(fwk$smp_data))
    stop("sub_sampling_var must be NULL, a column name of smp_data, or a numeric vector with one element per row of smp_data.", call. = FALSE)
  as.numeric(v)[ok]
}

## Out-of-fold predictions of the configuration `params`/`nrounds` taken from an xgb_tune(keep_oof = TRUE) result,
## after checking that the tune is the same model on the same data. Stops with the reason otherwise.
.xgb_tune_oof_hat <- function(tune, fwk, params, nrounds, transformation, rescale_weights, variance_y, domains) {
  o <- tune$oof
  if (is.null(o)) stop("pool_tune: the tune result has no out-of-fold predictions; rerun xgb_tune(keep_oof = TRUE).", call. = FALSE)
  bad <- function(...) stop("pool_tune does not match this xgb() call: ", ..., call. = FALSE)
  if (!identical(o$transformation, transformation)) bad("transformation '", o$transformation, "' vs '", transformation, "'.")
  if (!identical(as.character(o$cluster_var), as.character(domains))) bad("the tune's folds are grouped by '", o$cluster_var, "', not by domains = '", domains, "'.")
  if (!identical(isTRUE(o$rescale_weights), isTRUE(rescale_weights))) bad("rescale_weights differs.")
  if (!identical(o$variance_y, variance_y)) bad("variance_y differs.")
  if (!identical(as.character(o$covariates), as.character(fwk$covariates))) bad("covariates differ (or differ in order).")
  n <- nrow(fwk$smp_data)
  if (length(o$outcome) != n) bad("the tune has ", length(o$outcome), " rows, smp_data ", n, ".")
  if (!isTRUE(all.equal(as.numeric(o$outcome), as.numeric(fwk$Y_smp[[1]]), tolerance = 1e-12))) bad("the outcome differs.")
  if (!isTRUE(all.equal(as.numeric(o$weights), as.numeric(fwk$smp_weights_vec), tolerance = 1e-12))) bad("the survey weights differ.")
  if (!identical(as.character(o$domain), as.character(fwk$X_smp[[domains]]))) bad("the domain labels differ.")
  cand <- tune$candidates
  hp <- c("max_depth", "colsample_bytree", "colsample_bylevel", "colsample_bynode", "subsample", "min_child_weight",
          "eta", "gamma", "max_delta_step", "lambda", "alpha")
  hit <- rep(TRUE, nrow(cand))
  for (h in hp) hit <- hit & abs(cand[[h]] - params[[h]]) <= 1e-12 * pmax(1, abs(params[[h]]))
  hit <- hit & o$nround_used == nrounds
  if (!any(hit)) bad("no candidate of the tune has this configuration's hyperparameters and nrounds.")
  j <- which(hit)[1]
  o$pred[, j]
}

## The pools for each fit. `fits` are the point_estim_xgb results (one, or one per configuration); returns the
## fits with resid_sub_domains, resid_domains, wt_sub_domains and wt_domains replaced, plus a record.
.xgb_apply_pools <- function(fits, cfg_params, cfg_nrounds, fwk, X_smp_xgb, X_pop_xgb, transform_outcome, transformation,
                             rescale_weights, variance_y, domains, sub_domains, area_pool, sub_pool, area_deflate, sub_deflate,
                             area_sampling_var, sub_sampling_var, pool_folds, pool_seed, pool_fit_seed, pool_tune, seed) {
  cell <- fwk$smp_data[[sub_domains]]
  rows_pop <- match(cell, fwk$X_pop[[sub_domains]]); ok <- !is.na(rows_pop)
  if (any(vapply(fits, function(f) length(f$resid_total) != sum(ok), TRUE)))
    stop("out-of-fold pools: the sample rows do not map one to one onto census cells.", call. = FALSE)
  dom_smp <- as.character(fwk$X_smp[[domains]])
  dom_pop <- as.character(fwk$X_pop[[domains]][rows_pop])[ok]
  y_resp <- as.numeric(fwk$Y_smp[[1]]); w <- as.numeric(fwk$smp_weights_vec)
  y_t <- as.numeric(transform_outcome(y_resp)$y)
  het <- if (is.null(variance_y)) rep(1, length(y_t)) else as.numeric(fwk$smp_data[[variance_y]])^-0.5
  Xcell <- as.matrix(X_pop_xgb)[rows_pop, , drop = FALSE]
  K <- if (is.data.frame(pool_folds)) NA_integer_ else as.integer(pool_folds)
  fold_a <- fold_s <- NULL
  need_a <- area_pool == "out_of_fold" || sub_pool == "out_of_fold_domain"
  if (need_a)
    fold_a <- if (!is.null(pool_tune)) as.integer(pool_tune$oof$fold) else
      .xgb_pool_fold_domain(dom_smp, K, pool_seed, if (is.data.frame(pool_folds)) pool_folds else NULL)
  if (sub_pool == "out_of_fold") {
    if (is.na(K)) K <- 10L
    fold_s <- .xgb_pool_fold_sub(cell, dom_smp, K, pool_seed)
  }
  if (!is.null(pool_tune) && !isTRUE(all.equal(unname(as.matrix(X_smp_xgb)[ok, , drop = FALSE]), unname(Xcell[ok, , drop = FALSE]),
                                               tolerance = 1e-10)))
    stop("pool_tune: the tuning predictions are made at each sample row's own covariates, which differ from the covariates of its ",
         "census cell, so they are not the residuals xgb() defines. Use fresh folds (leave pool_tune NULL).", call. = FALSE)
  info <- vector("list", length(fits))
  for (k in seq_along(fits)) {
    f <- fits[[k]]
    pa <- list(resid = f$resid_domains, w = f$wt_domains, nm = f$domains)
    ps <- list(resid = f$resid_sub_domains, w = f$wt_sub_domains)
    if (need_a) {
      hat <- if (!is.null(pool_tune)) .xgb_tune_oof_hat(pool_tune, fwk, cfg_params[[k]], cfg_nrounds[k], transformation,
                                                         rescale_weights, variance_y, domains) else
        .xgb_oof_hat(fold_a, X_smp_xgb, y_t, w, dom_smp, het, Xcell, cfg_params[[k]], cfg_nrounds[k], rescale_weights, seed, pool_fit_seed)
      p_a <- .xgb_pools_from_resid((y_t - hat)[ok], w[ok], dom_pop)
      if (area_pool == "out_of_fold") pa <- list(resid = p_a$resid_domains, w = p_a$wt_domains, nm = p_a$names)
      if (sub_pool == "out_of_fold_domain") ps <- list(resid = p_a$resid_sub_domains, w = p_a$wt_sub_domains)
    }
    if (sub_pool == "out_of_fold") {
      hat <- .xgb_oof_hat(fold_s, X_smp_xgb, y_t, w, dom_smp, het, Xcell, cfg_params[[k]], cfg_nrounds[k], rescale_weights, seed, pool_fit_seed)
      p <- .xgb_pools_from_resid((y_t - hat)[ok], w[ok], dom_pop)
      ps <- list(resid = p$resid_sub_domains, w = p$wt_sub_domains)
    }
    fa <- fs <- NULL
    if (area_deflate) {
      v <- .xgb_resolve_area_var(area_sampling_var, fwk, ok, dom_pop, pa$nm, y_resp[ok], w[ok])
      est <- as.numeric(collapse:::fmean(y_resp[ok], g = dom_pop, w = w[ok], use.g.names = TRUE)[pa$nm])
      fa <- .xgb_deflate_pool(pa$resid, pa$w, .xgb_var_to_model_scale(v, est, transformation, transform_outcome), "area")
      pa$resid <- fa$resid
    }
    if (sub_deflate) {
      v <- .xgb_resolve_sub_var(sub_sampling_var, fwk, ok, variance_y)
      fs <- .xgb_deflate_pool(ps$resid, ps$w, .xgb_var_to_model_scale(v, y_resp[ok], transformation, transform_outcome), "sub")
      ps$resid <- fs$resid
    }
    f$resid_domains <- pa$resid; f$wt_domains <- pa$w; f$domains <- pa$nm
    f$resid_sub_domains <- ps$resid; f$wt_sub_domains <- ps$w
    fits[[k]] <- f
    info[[k]] <- list(resid_domains = pa$resid, wt_domains = pa$w, domains = pa$nm,
                      resid_sub_domains = ps$resid, wt_sub_domains = ps$w,
                      area_deflation = if (is.null(fa)) NULL else fa[c("factor", "msr", "msv", "n_missing")],
                      sub_deflation = if (is.null(fs)) NULL else fs[c("factor", "msr", "msv", "n_missing")])
  }
  list(fits = fits,
       record = list(area_pool = area_pool, sub_pool = sub_pool, area_deflate = area_deflate, sub_deflate = sub_deflate,
                     area_folds = if (is.null(fold_a)) NULL else data.frame(domain = unique(dom_smp), fold = fold_a[match(unique(dom_smp), dom_smp)]),
                     sub_folds = if (is.null(fold_s)) NULL else unique(data.frame(domain = dom_smp, sub_domain = as.character(cell), fold = fold_s)),
                     area_from_tune = !is.null(pool_tune), pool_seed = pool_seed, pool_fit_seed = pool_fit_seed, configs = info))
}
