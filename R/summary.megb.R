#' Summarize a Mixed Effects Gradient Boosting model
#'
#' Extracts information about the data, model fit, and random-effects
#' components from a \code{megb} object. The returned object is suitable for
#' printing with \code{print}.
#'
#' @param object an object of class \code{megb}, as returned by
#'   \code{\link{megb}}.
#' @param ... optional additional inputs (ignored).
#'
#' @return An object of class \code{c("summary.megb", "emdi")} with elements
#'   \code{call_megb}, \code{dom_info}, \code{smp_size}, \code{pop_size},
#'   \code{sizedom_smp_pop}, \code{coeff_determ}, \code{var_decomp}, and
#'   \code{megb_info}.
#'
#' @export
summary.megb <- function(object, ...) {

  call_megb <- object$out_call

  # ── Domain info ───────────────────────────────────────────────────────────────
  total_dom <- object$framework$domains_total
  in_dom    <- object$framework$domains_in
  oos_dom   <- object$framework$domains_out

  dom_info <- data.frame(in_dom, oos_dom, total_dom, row.names = "")
  colnames(dom_info) <- c("In-sample", "Out-of-sample", "Total")

  smp_size <- object$framework$N_smp
  pop_size <- object$framework$N_pop

  smp_size_dom <- summary(as.data.frame(object$framework$ni_smp)[, "Freq"])
  pop_size_dom <- summary(as.data.frame(object$framework$ni_pop)[, "Freq"])

  sizedom_smp_pop <- rbind(
    Sample_domains     = smp_size_dom,
    Population_domains = pop_size_dom
  )

  # ── Goodness-of-fit on sample ─────────────────────────────────────────────────
  # unit_pred_smp stores back-transformed sample predictions
  Y_actual <- object$framework$Y_smp
  Y_hat    <- object$framework$unit_pred_smp

  if (!is.null(Y_hat) && length(Y_hat) == length(Y_actual)) {
    r_squared <- tryCatch(cor(Y_actual, Y_hat)^2,    error = function(e) NA)
    mae       <- mean(abs(Y_actual - Y_hat))
    rank_cor  <- tryCatch(cor(Y_actual, Y_hat, method = "spearman"), error = function(e) NA)

    # Area-level metrics: weighted-mean of (Y_actual, Y_hat) within each domain,
    # then correlation / MAE across area means. Mirrors summary.xgb.
    dom_vec <- object$framework$smp_data[[object$framework$domains]]
    w_vec   <- object$framework$smp_weights_vec
    if (is.null(w_vec)) w_vec <- rep(1, length(Y_actual))
    yyhat   <- data.frame(Y = Y_actual, hat = Y_hat)
    area_means <- tryCatch(
      collapse::fmean(x = yyhat, g = dom_vec, w = w_vec),
      error = function(e) NULL
    )
    if (!is.null(area_means) && nrow(area_means) >= 2) {
      area_r_squared <- tryCatch(cor(area_means$Y, area_means$hat)^2,
                                 error = function(e) NA)
      area_mae       <- mean(abs(area_means$Y - area_means$hat))
      area_rank_cor  <- tryCatch(
        cor(area_means$Y, area_means$hat, method = "spearman"),
        error = function(e) NA
      )
    } else {
      area_r_squared <- NA; area_mae <- NA; area_rank_cor <- NA
    }
  } else {
    r_squared <- NA; mae <- NA; rank_cor <- NA
    area_r_squared <- NA; area_mae <- NA; area_rank_cor <- NA
  }

  coeff_det <- data.frame(
    R2            = r_squared,
    MAE           = mae,
    Rank_cor      = rank_cor,
    Area_R2       = area_r_squared,
    Area_MAE      = area_mae,
    Area_Rank_cor = area_rank_cor,
    row.names = ""
  )

  # ── Variance decomposition ────────────────────────────────────────────────────
  if (!is.null(Y_hat) && length(Y_hat) == length(Y_actual)) {
    domains_vec     <- object$framework$smp_data[[object$framework$domains]]
    resid           <- Y_actual - Y_hat
    domain_mean_res <- tapply(resid, domains_vec, mean)
    sigma2_v        <- var(domain_mean_res)
    resid_adj       <- resid - domain_mean_res[match(domains_vec,
                                                     names(domain_mean_res))]
    sigma2_e        <- var(resid_adj)
  } else {
    sigma2_v <- NA
    sigma2_e <- NA
  }

  var_decomp <- data.frame(
    Sigma2_v  = sigma2_v,
    Sigma2_e  = sigma2_e,
    row.names = ""
  )

  # ── Shrinkage factors ─────────────────────────────────────────────────────────
  ran_eff_sd   <- object$megb_model$ran_eff_sd
  error_sd_val <- object$megb_model$error_sd
  ni_smp_vec   <- as.numeric(as.data.frame(object$framework$ni_smp)[, "Freq"])

  if (!is.null(ran_eff_sd) && !is.na(ran_eff_sd) && ran_eff_sd > 0 &&
      !is.null(error_sd_val) && !is.na(error_sd_val) && error_sd_val > 0) {
    gamma_i  <- ran_eff_sd^2 / (ran_eff_sd^2 + error_sd_val^2 / ni_smp_vec)
  } else {
    gamma_i  <- rep(0, length(ni_smp_vec))
  }
  shrinkage <- summary(gamma_i)

  # ── MEGB model info ───────────────────────────────────────────────────────────
  megb_info <- data.frame(
    c(object$transformation,
      object$gbm_engine,
      if (!is.null(ran_eff_sd))   round(ran_eff_sd,   6) else NA,
      if (!is.null(error_sd_val)) round(error_sd_val, 6) else NA)
  )
  colnames(megb_info) <- NULL
  rownames(megb_info) <- c(
    "Transformation",
    "GB Engine",
    "Ran. Eff. SD",
    "Error SD"
  )

  sum_megb <- list(
    call_megb       = call_megb,
    dom_info        = dom_info,
    smp_size        = smp_size,
    pop_size        = pop_size,
    sizedom_smp_pop = sizedom_smp_pop,
    coeff_determ    = coeff_det,
    var_decomp      = var_decomp,
    shrinkage       = shrinkage,
    megb_info       = megb_info
  )

  class(sum_megb) <- c("summary.megb", "emdi")
  sum_megb
}

#' @export
print.summary.megb <- function(x, ...) {
  cat("________________________________________________________________\n")
  cat("Mixed Effects Gradient Boosting for Small Area Estimation\n")
  cat("________________________________________________________________\n")
  cat("Call:\n")
  print(x$call_megb)
  cat("\nDomains\n")
  cat("________________________________________________________________\n")
  print(x$dom_info)
  cat("\nTotals:\n")
  cat("Units in sample:", x$smp_size, "\n")
  if (!is.null(x$pop_size))
    cat("Units in population:", x$pop_size, "\n")
  cat("\n")
  print(x$sizedom_smp_pop)
  cat("\n")
  print(x$coeff_determ)
  cat("\nVariance components:\n")
  cat("________________________________________________________________\n")
  cat("Sigma2_v: area random-effect variance\n")
  cat("Sigma2_e: idiosyncratic error variance\n")
  print(x$var_decomp)
  cat("\nShrinkage factors (gamma_i):\n")
  cat("________________________________________________________________\n")
  sh <- matrix(x$shrinkage, nrow = 1,
               dimnames = list("Gamma", names(x$shrinkage)))
  print(sh)
  cat("\nMEGB model components:\n")
  cat("________________________________________________________________\n")
  print(x$megb_info)
  cat("\n")
}
