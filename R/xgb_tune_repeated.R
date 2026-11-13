#' Tune on repeated cross-validation, the shared estimation convention
#'
#' Draws \code{n_iter} configurations once at \code{draw_seed} from \code{bounds}
#' (default: the convention's bounds, see \code{xgb_tune}), then scores the same
#' configurations with \code{\link{xgb_tune}} on each fold seed in \code{fold_seeds}
#' (default 123, 124 and 125; ten domain-grouped folds). Single-seed selection is
#' mostly fold noise, so \code{\link{xgb_top_configs}} averages losses and paired
#' differences over the seeds. The result is the list of the \code{xgb_tune}
#' results, one per fold seed, all on the same candidates in the same order, which
#' \code{xgb_top_configs} takes as it is.
#'
#' @param fixed,smp_data,domains,... passed to \code{\link{xgb_tune}}. \code{grid},
#'   \code{search}, \code{seed}, \code{n_iter} and \code{bounds} are set here.
#' @param fold_seeds integer vector of fold seeds. Default \code{c(123, 124, 125)}.
#' @param draw_seed seed of the one draw of configurations. Default 20261010.
#' @param n_iter number of configurations. Default 150.
#' @param bounds named list of ranges as in \code{xgb_tune}; \code{NULL} takes the
#'   convention's.
#' @param fit_seed seed passed to xgboost for every fit; \code{NULL} uses \code{draw_seed}.
#' @return a list of \code{xgb_tune} results named \code{seed_<fold seed>}, with
#'   attributes \code{draw_seed}, \code{fold_seeds}, \code{fit_seed} and \code{bounds}.
#' @seealso \code{\link{xgb_tune}}, \code{\link{xgb_top_configs}}
#' @export
xgb_tune_repeated <- function(fixed, smp_data, domains, fold_seeds = c(123, 124, 125), draw_seed = 20261010,
                              n_iter = 150, bounds = NULL, fit_seed = NULL, ...) {
  if (is.null(bounds)) bounds <- .XGB_DEFAULT_BOUNDS
  stopifnot(is.list(bounds), length(fold_seeds) >= 1, !anyDuplicated(fold_seeds), n_iter >= 1)
  if (is.null(fit_seed)) fit_seed <- draw_seed
  had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  if (had) old <- get(".Random.seed", envir = globalenv())
  on.exit(if (had) assign(".Random.seed", old, envir = globalenv()) else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) rm(".Random.seed", envir = globalenv()), add = TRUE)
  set.seed(draw_seed)
  cand <- .xgb_sample_bounds(bounds, n_iter)
  out <- lapply(fold_seeds, function(s)
    xgb_tune(fixed = fixed, smp_data = smp_data, domains = domains, grid = cand, seed = s, fit_seed = fit_seed,
             keep_candidates = TRUE, ...))
  names(out) <- paste0("seed_", fold_seeds)
  attr(out, "draw_seed") <- draw_seed; attr(out, "fold_seeds") <- fold_seeds
  attr(out, "fit_seed") <- fit_seed; attr(out, "bounds") <- bounds
  out
}
