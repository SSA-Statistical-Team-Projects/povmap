# Internal MSE bootstrap helpers for megb_em
#' @importFrom dplyr left_join group_by summarise
#' @importFrom lme4 ranef

# Block-sample errors by domain. When weights are supplied AND weightedBS is
# TRUE, residuals are drawn with probability proportional to weights (matching
# xgb's weightedBS semantics). When weights are NULL or weightedBS is FALSE,
# residuals are drawn with equal probability (historical behaviour).
block_sample <- function(domains, in_samp, smp_data, dom_name, pop_data, gb_res,
                          weights = NULL, weightedBS = FALSE) {
  block_err <- vector(mode = "list", length = length(domains))
  use_w     <- !is.null(weights) && isTRUE(weightedBS)

  for (idd in which(in_samp)) {
    in_dom <- which(smp_data[[dom_name]] == domains[idd])
    x_vals <- gb_res[in_dom]
    p_vals <- if (use_w) weights[in_dom] / sum(weights[in_dom]) else NULL
    block_err[[idd]] <- sample(
      x_vals,
      size    = sum(pop_data[dom_name] == domains[idd]),
      replace = TRUE,
      prob    = p_vals
    )
  }

  if (sum(in_samp) != length(domains)) {
    # OOS domains draw from the full residual pool, weighted across the entire
    # in-sample if weights are in use.
    pool_prob <- if (use_w) weights / sum(weights) else NULL
    for (idd in which(!in_samp)) {
      block_err[[idd]] <- sample(
        gb_res,
        size    = sum(pop_data[dom_name] == domains[idd]),
        replace = TRUE,
        prob    = pool_prob
      )
    }
  }

  unlist(block_err)
}

# Stratified sample selection for bootstrap
sample_select <- function(pop, smp, dom_name) {
  smpSizes <- table(smp[dom_name])
  smpSizes <- data.frame(
    smpidD = as.character(names(smpSizes)),
    n_smp  = as.numeric(smpSizes),
    stringsAsFactors = FALSE
  )

  smpSizes <- dplyr::left_join(
    data.frame(idD = as.character(unique(pop[[dom_name]]))),
    smpSizes,
    by = c("idD" = "smpidD")
  )
  smpSizes$n_smp[is.na(smpSizes$n_smp)] <- 0

  splitPop <- split(pop, pop[[dom_name]])

  stratSamp <- function(dfList, ns) {
    do.call(rbind, mapply(dfList, ns, FUN = function(df, n) {
      sel <- base::sample(seq_len(nrow(df)), n, replace = n > nrow(df))
      df[sel, ]
    }, SIMPLIFY = FALSE))
  }

  stratSamp(splitPop, smpSizes$n_smp)
}

# Compute random components for bootstrap resampling
#
# DESIGN NOTE — why we do NOT domain-demean the residuals:
#
# The old code computed ran_effs as domain means of (Y - unit_pred_smp) and
# gb_res as the within-domain deviations. Because unit_pred_smp already
# includes the BLUP û_d, these domain means are ≈ 0 for all domains, so
# ran_effs was just rescaled noise (not actual BLUPs). Worse, for singleton
# wards (n_d = 1), the within-domain deviation gb_eij = 0 exactly (the single
# residual minus itself). With many singleton wards this diluted sd(gb_eij),
# causing the rescaling to (gb_eij / sd(gb_eij)) * σ_e to badly over-inflate
# the multi-obs residuals — which is why σ_e_boot ≈ 0.24 >> σ_e = 0.157 in
# the bootstrap, and σ_u_boot was correspondingly deflated.
#
# Correct approach (Prasad-Rao / Hall-Maiti parametric bootstrap):
#   ran_effs = actual û_d BLUPs from the fitted lme4 (via lme4::ranef),
#              rescaled slightly to σ_u to correct for finite-sample shrinkage.
#   gb_res   = Y - unit_pred_smp = Y - gb_smp - û_d ≈ e_dk
#              (pure idiosyncratic residuals, not domain-demeaned).
#
ran_comp <- function(unit_pred_smp, smp_data, Y, dom_name, error_sd,
                     cov_names = cov_names, model, weights = NULL) {
  # weights: numeric vector aligned to rows of smp_data. When supplied,
  # gb_res rescaling uses weighted SD and weighted-mean centering so the
  # residual pool is properly normalised under the survey design.
  smp_data_tmp <- smp_data
  residuals_gb <- Y - unit_pred_smp   # ≈ e_dk (idiosyncratic residuals)

  # Helpers: weighted SD / mean fall back to unweighted when weights are NULL.
  wmean <- function(x, w) {
    if (is.null(w)) mean(x) else stats::weighted.mean(x, w)
  }
  wsd <- function(x, w) {
    if (is.null(w)) sd(x)
    else {
      m   <- stats::weighted.mean(x, w)
      sqrt(sum(w * (x - m)^2) / sum(w))
    }
  }

  if (is.null(model$ran_eff_sd) || model$ran_eff_sd < 1e-8) {
    # Singular LMM: all random effects are zero; use raw residuals for e*.
    gb_res_sd <- wsd(residuals_gb, weights)
    gb_res    <- if (gb_res_sd > 1e-10) (residuals_gb / gb_res_sd) * error_sd
                 else residuals_gb
    gb_res    <- gb_res - wmean(gb_res, weights)
    ran_effs  <- rep(0, length(unique(smp_data[[dom_name]])))
  } else {
    # ── ran_effs: extract actual per-domain BLUPs û_d from the fitted lme4 ──
    blup_df  <- as.data.frame(lme4::ranef(model$effect_model)[[dom_name]])
    blup_vec <- setNames(blup_df[, "(Intercept)"], rownames(blup_df))

    insamp_doms <- as.character(unique(smp_data[[dom_name]]))
    ran_effs    <- as.numeric(blup_vec[insamp_doms])
    ran_effs[is.na(ran_effs)] <- 0L

    # Rescale to model$ran_eff_sd. ran_effs is per-domain (length = D, not n),
    # so weights at the area level (sum of within-area sample weights) are the
    # right reference if we want a weighted rescaling. Without that, fall back
    # to unweighted SD across domains (the historical behaviour).
    area_w <- if (!is.null(weights)) {
      as.numeric(tapply(weights, smp_data[[dom_name]], sum)[insamp_doms])
    } else NULL
    re_sd <- wsd(ran_effs, area_w)
    if (re_sd > 1e-10) ran_effs <- (ran_effs / re_sd) * model$ran_eff_sd
    ran_effs <- ran_effs - wmean(ran_effs, area_w)

    # ── gb_res: idiosyncratic residuals, NOT domain-demeaned ─────────────────
    gb_res_sd <- wsd(residuals_gb, weights)
    gb_res    <- if (gb_res_sd > 1e-10) (residuals_gb / gb_res_sd) * error_sd
                 else residuals_gb
    gb_res    <- gb_res - wmean(gb_res, weights)
  }

  smp_data_tmp$gb_res <- gb_res
  list(gb_res = gb_res, ran_effs = ran_effs, smp_data = smp_data_tmp,
       area_w = if (!is.null(weights))
         as.numeric(tapply(weights, smp_data[[dom_name]], sum)[as.character(unique(smp_data[[dom_name]]))])
         else NULL)
}
