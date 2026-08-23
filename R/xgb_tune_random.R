## Helpers for xgb_tune(search = "random").
##
## Random search samples from CONTINUOUS RANGES rather than from a coarse
## candidate pool, because sampling a pool just reproduces the limitation of the
## grid it replaces. Parameters spanning orders of magnitude are drawn
## log-uniformly, since the difference between eta 0.01 and 0.05 matters far more
## than between 0.25 and 0.29 and a uniform draw spends most of its budget at the
## insensitive end.

## Default search ranges. Only these six vary unless the caller names more in
## `bounds`; a user exploring depth should not silently get gamma and alpha
## moving as well.
.XGB_DEFAULT_BOUNDS <- list(
  eta              = c(0.01, 0.3),
  max_depth        = c(2, 10),
  min_child_weight = c(1, 50),
  subsample        = c(0.5, 1.0),
  colsample_bytree = c(0.3, 1.0),
  lambda           = c(0.1, 10),
  ## nround is searched rather than fixed by early stopping: unit-level
  ## stopping underfits a domain-level objective (measured: patience 20
  ## halted at 34 rounds when the optimum was near 1258, halving the r2).
  nround           = c(50, 1500)
)

## Default sampling scale per parameter.
##   "log"     - spans orders of magnitude
##   "int"     - integer valued
##   "zeroinf" - non-negative with a point mass at zero; log is impossible at 0
##               and a plain uniform spends the budget on values nobody wants
##   "linear"  - everything else
.XGB_DEFAULT_SCALE <- list(
  eta               = "log",
  min_child_weight  = "log",
  lambda            = "log",
  max_depth         = "int",
  alpha             = "zeroinf",
  gamma             = "linear",
  subsample         = "linear",
  colsample_bytree  = "linear",
  colsample_bylevel = "linear",
  colsample_bynode  = "linear",
  max_delta_step    = "linear",
  nround            = "int"
)

#' Normalise one bounds entry
#'
#' Accepts either the short form \code{c(lo, hi)}, which inherits the
#' parameter's default scale, or the explicit form
#' \code{list(range = c(lo, hi), scale = "log")}. The short form deliberately
#' avoids \code{c(0.01, 0.3, "log")}, which would coerce the numbers to
#' character.
#' @keywords internal
.xgb_norm_bound <- function(param, b) {
  default_scale <- .XGB_DEFAULT_SCALE[[param]]
  if (is.null(default_scale)) default_scale <- "linear"
  if (is.list(b)) {
    rng <- b$range
    scl <- if (is.null(b$scale)) default_scale else b$scale
    p0  <- if (is.null(b$p0)) 0.5 else b$p0
  } else {
    rng <- b
    scl <- default_scale
    p0  <- 0.5
  }
  if (length(rng) != 2L || !is.numeric(rng) || any(!is.finite(rng)))
    stop("bounds for '", param, "' must be a numeric c(lo, hi).", call. = FALSE)
  if (rng[2] < rng[1])
    stop("bounds for '", param, "': upper limit is below the lower limit.", call. = FALSE)
  if (scl == "log" && rng[1] <= 0)
    stop("bounds for '", param, "' use a log scale, so the lower limit must be > 0.",
         call. = FALSE)
  list(range = rng, scale = scl, p0 = p0)
}

#' Draw n configurations from the bounds
#' @keywords internal
.xgb_sample_bounds <- function(bounds, n) {
  cols <- lapply(names(bounds), function(p) {
    s  <- .xgb_norm_bound(p, bounds[[p]])
    lo <- s$range[1]; hi <- s$range[2]
    if (s$scale == "log") {
      exp(stats::runif(n, log(lo), log(hi)))
    } else if (s$scale == "int") {
      if (hi == lo) rep(as.integer(lo), n)
      else sample(seq.int(as.integer(lo), as.integer(hi)), n, replace = TRUE)
    } else if (s$scale == "zeroinf") {
      ifelse(stats::runif(n) < s$p0, 0, stats::runif(n, lo, hi))
    } else {
      stats::runif(n, lo, hi)
    }
  })
  names(cols) <- names(bounds)
  as.data.frame(cols, stringsAsFactors = FALSE)
}

#' Where does the selected value sit inside its range?
#'
#' Returns the position as a fraction of the range together with a boundary
#' flag. An optimum sitting on an edge means the range was too narrow, and the
#' caller should widen it rather than have the search silently truncate.
#' @keywords internal
.xgb_bounds_position <- function(bounds, chosen, tol = 0.05) {
  rows <- lapply(names(bounds), function(p) {
    s <- .xgb_norm_bound(p, bounds[[p]])
    lo <- s$range[1]; hi <- s$range[2]; v <- chosen[[p]]
    pos <- if (hi > lo) {
      if (s$scale == "log") (log(v) - log(lo)) / (log(hi) - log(lo))
      else (v - lo) / (hi - lo)
    } else 0.5
    data.frame(parameter = p, lower = lo, upper = hi, scale = s$scale,
               selected = v, position = pos,
               at_boundary = isTRUE(pos <= tol || pos >= 1 - tol),
               stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

## ---------------------------------------------------------------------------
## Shared scoring
##
## Grid, random and Bayesian search all score a configuration the same way:
## grouped CV holding out whole clusters, predictions aggregated to weighted
## domain means, mean squared error at the DOMAIN level. Keeping that in one
## place is the point -- reimplementing it per optimiser is how the metric
## drifts between search modes.
## ---------------------------------------------------------------------------

#' Score every row of a tuning grid by grouped cross-validation
#'
#' Assumes a foreach parallel backend is already registered by the caller.
#' Work is flattened over (fold, configuration) pairs so a single-row grid, as
#' produced by each Bayesian optimisation step, still uses the whole backend.
#'
#' @return list with `OPT` and `ITER`, both configurations x folds, and
#'   `domain_labels`, the held-out domain means per fold
#' @keywords internal
.xgb_score_configs <- function(tunegrid, X_final, X_smp, Y_smp, smp_weights,
                               cluster_col, cluster, domains, folds,
                               early_stopping_rounds = NULL, nrounds_max = 2000,
                               seed = NULL, dots = list()) {
  nrow_g <- nrow(tunegrid)
  es_on  <- !is.null(early_stopping_rounds)
  OPT    <- matrix(NA_real_, nrow = nrow_g, ncol = folds)
  ITER   <- matrix(NA_real_, nrow = nrow_g, ncol = folds)
  domain_labels_list <- vector("list", folds)

  ## per-fold quantities that do not depend on the configuration
  fold_info <- vector("list", folds)
  for (fold in seq_len(folds)) {
    hold <- cluster_col$fold == fold
    fold_df <- data.frame(labels = Y_smp[hold, ],
                          domains = X_smp[hold, ][[paste0(domains)]],
                          wts = smp_weights[hold])
    domain_labels_list[[fold]] <- sapply(split(fold_df, fold_df$domains),
                                         function(g) stats::weighted.mean(g$labels, w = g$wts))
    train_rows <- !hold
    if (es_on) {
      ## Inner early-stopping group, carved from the TRAINING clusters so the
      ## held-out fold stays purely for scoring.
      tr_clusters  <- unique(cluster_col[[cluster]][train_rows])
      val_clusters <- sample(tr_clusters, max(1L, round(0.1 * length(tr_clusters))))
      es_val   <- train_rows & (cluster_col[[cluster]] %in% val_clusters)
      es_train <- train_rows & !es_val
    } else {
      es_val <- NULL; es_train <- train_rows
    }
    fold_info[[fold]] <- list(hold = hold, es_train = es_train, es_val = es_val)
  }

  tasks <- expand.grid(row = seq_len(nrow_g), fold = seq_len(folds))
  res <- foreach::foreach(t = seq_len(nrow(tasks)), .combine = rbind,
                          .packages = "xgboost") %dopar% {
    row <- tasks$row[t]; fold <- tasks$fold[t]
    fi <- fold_info[[fold]]
    params <- c(list(
      max_depth         = tunegrid$max_depth[row],
      colsample_bytree  = tunegrid$colsample_bytree[row],
      colsample_bylevel = tunegrid$colsample_bylevel[row],
      subsample         = tunegrid$subsample[row],
      min_child_weight  = tunegrid$min_child_weight[row],
      eta               = tunegrid$eta[row],
      gamma             = tunegrid$gamma[row],
      max_delta_step    = tunegrid$max_delta_step[row],
      lambda            = tunegrid$lambda[row],
      alpha             = tunegrid$alpha[row]), dots)
    if (!is.null(seed)) params$seed <- seed

    dtrain <- xgboost::xgb.DMatrix(data = data.matrix(X_final[fi$es_train, ]),
                                   label = Y_smp[fi$es_train, ],
                                   weight = as.matrix(smp_weights)[fi$es_train, ])
    if (es_on) {
      dvalid <- xgboost::xgb.DMatrix(data = data.matrix(X_final[fi$es_val, ]),
                                     label = Y_smp[fi$es_val, ],
                                     weight = as.matrix(smp_weights)[fi$es_val, ])
      fit <- xgboost::xgb.train(data = dtrain, params = params, nrounds = nrounds_max,
                                evals = list(valid = dvalid),
                                early_stopping_rounds = early_stopping_rounds, verbose = 0)
      bi <- suppressWarnings(as.integer(xgboost::xgb.attr(fit, "best_iteration")))
      best_iter <- if (length(bi) == 1L && !is.na(bi)) bi + 1L else nrounds_max
    } else {
      fit <- xgboost::xgb.train(data = dtrain, params = params,
                                nrounds = tunegrid$nround[row], verbose = 0)
      best_iter <- tunegrid$nround[row]
    }

    hat <- predict(fit, data.matrix(X_final[fi$hold, ]))
    dh  <- data.frame(hat = hat,
                      domains = X_smp[fi$hold, ][[paste0(domains)]],
                      labels = Y_smp[fi$hold, ],
                      wts = smp_weights[fi$hold])
    g   <- split(dh, dh$domains)
    mh  <- sapply(g, function(z) stats::weighted.mean(z$hat, w = z$wts))
    ml  <- sapply(g, function(z) stats::weighted.mean(z$labels, w = z$wts))
    c(row = row, fold = fold, mse = mean((ml - mh)^2), best_iter = best_iter)
  }
  res <- matrix(as.numeric(res), ncol = 4)
  for (k in seq_len(nrow(res))) {
    OPT[res[k, 1], res[k, 2]]  <- res[k, 3]
    ITER[res[k, 1], res[k, 2]] <- res[k, 4]
  }
  list(OPT = OPT, ITER = ITER, domain_labels = domain_labels_list)
}
