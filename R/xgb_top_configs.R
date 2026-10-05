#' Select the configurations to average from one or more xgb_tune results
#'
#' Configuration averaging fits several near-equivalent tuned configurations
#' and averages their domain estimates (see the \code{configs} argument of
#' \code{\link{xgb}}). This function picks the set from the cross-validation
#' results of \code{\link{xgb_tune}}.
#'
#' @param tune the result of \code{\link{xgb_tune}} run with
#'   \code{keep_candidates = TRUE}, or a list of such results on the same grid
#'   but with different fold splits (for example, different \code{seed}s). With
#'   several splits each configuration is scored by its mean cross-validation
#'   error over all of them, which makes the set much less sensitive to the
#'   fold assignment.
#' @param rule how the set is chosen around the best configuration (lowest
#'   score). \code{"one_se_paired"} (default): every configuration whose mean
#'   fold-wise difference from the best is at most the standard error of that
#'   difference. Comparing fold by fold removes the variation the
#'   configurations share. With several splits the difference is averaged over
#'   splits and its standard error is the typical single-split one (it is not
#'   divided further by the number of splits, because the splits reuse the same
#'   sample). \code{"one_se"}: every configuration whose score is at most the
#'   best score plus the best configuration's own standard error.
#'   \code{"k_best"}: the \code{k} configurations with the lowest scores.
#' @param k number of configurations for \code{rule = "k_best"}.
#' @param floor,cap the set is filled up to \code{floor} configurations with the
#'   next best by score, and cut to the \code{cap} best by score. Defaults 3 and 8.
#' @return a data.frame with one row per selected configuration, best first:
#'   the twelve hyperparameters under the names \code{\link{xgb}} uses
#'   (\code{nrounds}, \code{max_depth}, \code{colsample_bytree},
#'   \code{colsample_bylevel}, \code{colsample_bynode}, \code{subsample},
#'   \code{min_child_weight}, \code{eta}, \code{gamma}, \code{max_delta_step},
#'   \code{lambda}, \code{alpha}), \code{weight} (equal weights summing to 1),
#'   \code{cv_score} (mean cross-validation error) and \code{rank} (rank of the
#'   score among all candidates). It can be passed to \code{xgb(configs = )}
#'   as it is.
#' @seealso \code{\link{xgb_tune}}, \code{\link{xgb}}
#' @export
xgb_top_configs <- function(tune, rule = c("one_se_paired", "one_se", "k_best"),
                            k = 5, floor = 3, cap = 8) {
  rule <- match.arg(rule)
  tunes <- if (!is.null(tune$candidates) || !is.null(tune$cv_by_fold)) list(tune) else tune
  if (!is.list(tunes) || !length(tunes))
    stop("tune must be an xgb_tune result or a list of xgb_tune results.", call. = FALSE)
  hp <- .xgb_hp_names
  for (t in tunes) {
    if (is.null(t$candidates) || is.null(t$cv_by_fold))
      stop("each xgb_tune result needs candidates and cv_by_fold; rerun xgb_tune with keep_candidates = TRUE.",
           call. = FALSE)
    if (nrow(t$candidates) != nrow(t$cv_by_fold))
      stop("candidates and cv_by_fold have different numbers of rows.", call. = FALSE)
  }
  grid <- tunes[[1]]$candidates
  if (!all(hp %in% names(grid))) stop("candidates lack hyperparameter columns: ",
                                      paste(setdiff(hp, names(grid)), collapse = ", "), call. = FALSE)
  for (t in tunes[-1])
    if (!isTRUE(all.equal(t$candidates[, hp], grid[, hp], check.attributes = FALSE)))
      stop("all xgb_tune results must be on the same grid, in the same order.", call. = FALSE)
  stopifnot(k >= 1, floor >= 1, cap >= floor)

  cvs <- lapply(tunes, function(t) as.matrix(t$cv_by_fold))
  score <- Reduce(`+`, lapply(cvs, rowMeans)) / length(cvs)
  ok <- is.finite(score)
  if (!any(ok)) stop("no configuration has a finite cross-validation score.", call. = FALSE)
  ord <- order(score, na.last = NA)          # finite scores, best first (ties keep grid order)
  b <- ord[1]
  inset <- switch(rule,
    k_best = ord[seq_len(min(k, length(ord)))],
    one_se = {
      se_b <- sqrt(mean(vapply(cvs, function(cf) stats::var(cf[b, ]) / ncol(cf), 1)))
      ord[score[ord] <= score[b] + se_b]
    },
    one_se_paired = {
      d  <- Reduce(`+`, lapply(cvs, function(cf) rowMeans(sweep(cf, 2, cf[b, ])))) / length(cvs)
      se <- sqrt(Reduce(`+`, lapply(cvs, function(cf) apply(sweep(cf, 2, cf[b, ]), 1, stats::var) / ncol(cf))) / length(cvs))
      keep <- ok & d <= se; keep[b] <- TRUE
      ord[keep[ord]]
    })
  if (length(inset) < floor) inset <- c(inset, setdiff(ord, inset)[seq_len(min(floor, length(ord)) - length(inset))])
  inset <- inset[order(score[inset])]
  inset <- inset[seq_len(min(cap, length(inset)))]

  out <- grid[inset, hp, drop = FALSE]
  names(out)[names(out) == "nround"] <- "nrounds"
  out$weight <- 1 / nrow(out)
  out$cv_score <- score[inset]
  out$rank <- match(inset, ord)
  rownames(out) <- NULL
  out
}

## xgb()'s twelve hyperparameter arguments and their defaults, read from its formals
.xgb_xgb_defaults <- function() {
  f <- formals(xgb)
  hp <- c("nrounds", "max_depth", "colsample_bytree", "colsample_bylevel", "colsample_bynode",
          "subsample", "min_child_weight", "eta", "gamma", "max_delta_step", "lambda", "alpha")
  stats::setNames(lapply(hp, function(h) eval(f[[h]])), hp)
}

## the twelve hyperparameters, as xgb_tune names them
.xgb_hp_names <- c("nround", "max_depth", "colsample_bytree", "colsample_bylevel", "colsample_bynode",
                   "subsample", "min_child_weight", "eta", "gamma", "max_delta_step", "lambda", "alpha")

## Validate and complete xgb(configs = ): one row per configuration with the twelve hyperparameters
## under xgb's names (nrounds or nround accepted), missing ones filled from `defaults`, and a weight
## (equal if absent) normalised to sum to 1. Columns cv_score, rank, mse and config are ignored.
.xgb_check_configs <- function(configs, defaults) {
  if (!is.data.frame(configs) || nrow(configs) < 1)
    stop("configs must be a data.frame with one row per configuration.", call. = FALSE)
  if ("nround" %in% names(configs) && !"nrounds" %in% names(configs)) names(configs)[names(configs) == "nround"] <- "nrounds"
  hp <- names(defaults)
  extra <- setdiff(names(configs), c(hp, "weight", "cv_score", "rank", "mse", "config"))
  if (length(extra)) stop("configs has unknown columns: ", paste(extra, collapse = ", "), call. = FALSE)
  out <- data.frame(row.names = seq_len(nrow(configs)))
  for (h in hp) {
    v <- if (h %in% names(configs)) configs[[h]] else rep(defaults[[h]], nrow(configs))
    if (!is.numeric(v) || anyNA(v)) stop("configs$", h, " must be numeric without NA.", call. = FALSE)
    out[[h]] <- as.numeric(v)
  }
  w <- if ("weight" %in% names(configs)) as.numeric(configs$weight) else rep(1, nrow(configs))
  if (anyNA(w) || any(w < 0) || sum(w) <= 0) stop("configs$weight must be non-negative with a positive sum.", call. = FALSE)
  out$weight <- w / sum(w)
  out
}
