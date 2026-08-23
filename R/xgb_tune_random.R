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
  lambda           = c(0.1, 10)
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
