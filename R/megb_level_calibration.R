#' Cross-validated level calibration of the MSE of unsampled areas (arm D)
#'
#' Implements arm D of MEGB_BOOT_FIX_PREREG_v8_final.md on a \code{megb()} fit made with \code{boot_synthetic = TRUE}.
#' The out-of-sample error of the cross-fitted synthetic prediction is measured on the sampled areas against a direct
#' estimate, \eqn{\hat M = \sum_i a_i[(y^{dir}_i - p^{cf}_i)^2 - v^{dir}_i] / \sum_i a_i}; the bootstrap's estimate of the
#' same quantity is \eqn{\hat B = \sum_i a_i \bar{(p^{*}_{ib} - \tau^{*}_{ib})^2} / \sum_i a_i} (option A of the
#' preregistration, from the replicates' synthetic predictions on the reporting scale); the calibration factor is
#' \eqn{\hat\kappa = \hat M / \hat B}, two-sided. The MSE of every unsampled area is multiplied by \eqn{\hat\kappa}; the
#' sampled areas' estimates, MSE and intervals, and all point estimates, are returned unchanged.
#'
#' \eqn{\hat\kappa} is computed from this fit alone (one draw, one number); it is never pooled across fits. It uses the
#' unbenchmarked synthetic predictions and is applied to whichever MSE \code{which} selects (benchmarked or not), which is
#' an assumption.
#'
#' @param fit a \code{megb} fit made with \code{mse = TRUE}, \code{bootstrap_refit = "leaves_only"} and \code{boot_synthetic = TRUE}.
#' @param direct data frame with columns \code{Domain}, \code{ydir} (direct estimate on the reporting scale) and
#'   \code{vdir} (its estimated variance). Areas that are not sampled, or have a missing \code{ydir} or \code{vdir}, are left
#'   out of \eqn{\hat M} and \eqn{\hat B} (their number is returned).
#' @param bounds length-two vector, the clipping bounds of \eqn{\hat\kappa}; default \code{c(0.1, 10)}.
#' @param weights \code{"equal"} (default, \eqn{a_i = 1}) or \code{"precision"} (\eqn{a_i = 1 / (v^{dir}_i + \hat B_i)}, a sensitivity check).
#' @param which \code{"bench"} (default when the fit is benchmarked) to scale and report the benchmarked MSE and interval, or \code{"unbench"}.
#' @param mse_type \code{"mse"} (default): the interval is \eqn{\mathrm{Mean} \pm z\sqrt{MSE}} as in \code{megb(mse_type = "mse")}; \code{"var"}: the legacy interval
#'   half-widths are rescaled by \eqn{\sqrt{\hat\kappa}} around the point estimate (an approximation where the legacy interval was clamped).
#' @param clamp optional length-two vector \code{c(lower, upper)} applied to the intervals (e.g. \code{c(0, 1)} for a proportion).
#' @param conf_level confidence level of the rebuilt interval (default 0.95; \code{mse_type = "mse"} only).
#' @return a list with \code{kappa} (after clipping and fallback), \code{kappa_raw}, \code{M_hat}, \code{B_hat}, \code{clipped} (\code{"lower"}, \code{"upper"} or \code{"none"}),
#'   \code{fallback}, \code{n_areas}, \code{n_left_out}, and \code{table}: one row per area with \code{Domain}, \code{sampled}, \code{Mean}, \code{MSE_legacy}, \code{MSE_D},
#'   \code{Lower_legacy}, \code{Upper_legacy}, \code{Lower_D}, \code{Upper_D}.
#' @export
megb_level_calibration <- function(fit, direct, bounds = c(0.1, 10), weights = c("equal", "precision"),
                                   which = c("bench", "unbench"), mse_type = c("mse", "var"), clamp = NULL, conf_level = 0.95) {
  weights  <- match.arg(weights); which <- match.arg(which); mse_type <- match.arg(mse_type)
  syn <- fit$boot_diag$synthetic
  if (is.null(syn)) stop("megb_level_calibration needs a fit made with mse = TRUE, bootstrap_refit = \"leaves_only\" and boot_synthetic = TRUE.")
  if (length(bounds) != 2L || !all(is.finite(bounds)) || bounds[1] <= 0 || bounds[1] >= bounds[2]) stop("bounds must be two positive numbers, lower < upper.")
  if (!all(c("Domain", "ydir", "vdir") %in% names(direct))) stop("direct needs the columns Domain, ydir and vdir.")
  dom <- as.character(syn$domains)
  d   <- data.frame(Domain = as.character(direct$Domain), ydir = as.numeric(direct$ydir), vdir = as.numeric(direct$vdir), stringsAsFactors = FALSE)
  i   <- match(d$Domain, dom)
  if (anyNA(i)) stop("direct contains domains that are not in the fit: ", paste(utils::head(d$Domain[is.na(i)], 5), collapse = ", "))
  use <- syn$sampled[i] & is.finite(d$ydir) & is.finite(d$vdir) & is.finite(syn$point[i]) & is.finite(syn$mse_synth[i])
  n_left_out <- sum(syn$sampled) - sum(use)   # every sampled area not used: missing from direct, or without a finite ydir, vdir or synthetic prediction
  if (!any(use)) stop("no sampled area with a finite direct estimate, variance and synthetic prediction.")
  e2   <- (d$ydir[use] - syn$point[i][use])^2 - d$vdir[use]
  Bi   <- syn$mse_synth[i][use]
  a    <- if (weights == "equal") rep(1, sum(use)) else 1 / (d$vdir[use] + Bi)
  M_hat <- sum(a * e2) / sum(a); B_hat <- sum(a * Bi) / sum(a)
  fallback <- !(is.finite(M_hat) && is.finite(B_hat) && M_hat > 0 && B_hat > 0)
  kappa_raw <- if (fallback) NA_real_ else M_hat / B_hat
  kappa <- if (fallback) 1 else min(max(kappa_raw, bounds[1]), bounds[2])
  clipped <- if (fallback) "none" else if (kappa_raw < bounds[1]) "lower" else if (kappa_raw > bounds[2]) "upper" else "none"
  # the legacy MSE and interval to scale
  vv <- fit$var; ci <- fit$CI; ind <- fit$ind
  if (which == "bench") {
    if (is.null(vv$Var_bench) || is.null(ind$Mean_bench)) stop("which = \"bench\": the fit has no benchmarked estimates; use which = \"unbench\".")
    mse0 <- vv$Var_bench[match(dom, as.character(vv$Domain))]; ctr <- ind$Mean_bench[match(dom, as.character(ind$Domain))]
    lo0 <- ci$Lower_bench[match(dom, as.character(ci$Domain))]; hi0 <- ci$Upper_bench[match(dom, as.character(ci$Domain))]
  } else {
    mse0 <- vv$Mean[match(dom, as.character(vv$Domain))]; ctr <- ind$Mean[match(dom, as.character(ind$Domain))]
    lo0 <- ci$Lower[match(dom, as.character(ci$Domain))]; hi0 <- ci$Upper[match(dom, as.character(ci$Domain))]
  }
  uns <- !syn$sampled
  mseD <- ifelse(uns, kappa * mse0, mse0)
  if (mse_type == "mse") {
    z <- stats::qnorm(1 - (1 - conf_level) / 2)
    loD <- ifelse(uns, ctr - z * sqrt(pmax(mseD, 0)), lo0); hiD <- ifelse(uns, ctr + z * sqrt(pmax(mseD, 0)), hi0)
  } else {
    loD <- if (kappa == 1) lo0 else ifelse(uns, ctr - sqrt(kappa) * (ctr - lo0), lo0)
    hiD <- if (kappa == 1) hi0 else ifelse(uns, ctr + sqrt(kappa) * (hi0 - ctr), hi0)
  }
  if (!is.null(clamp)) { loD <- ifelse(uns, pmax(loD, clamp[1]), loD); hiD <- ifelse(uns, pmin(hiD, clamp[2]), hiD) }
  list(kappa = kappa, kappa_raw = kappa_raw, M_hat = M_hat, B_hat = B_hat, clipped = clipped, fallback = fallback,
       n_areas = sum(use), n_left_out = n_left_out, weights = weights, which = which,
       table = data.frame(Domain = dom, sampled = syn$sampled, Mean = ctr, MSE_legacy = mse0, MSE_D = mseD,
                          Lower_legacy = lo0, Upper_legacy = hi0, Lower_D = loD, Upper_D = hiD, stringsAsFactors = FALSE))
}
