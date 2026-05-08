# Internal MSE bootstrap helpers for megb_em
#' @importFrom dplyr left_join group_by summarise
#' @importFrom lme4 ranef

# Block-sample errors by domain
block_sample <- function(domains, in_samp, smp_data, dom_name, pop_data, gb_res) {
  block_err <- vector(mode = "list", length = length(domains))

  for (idd in which(in_samp)) {
    block_err[[idd]] <- sample(
      gb_res[smp_data[dom_name] == domains[idd]],
      size    = sum(pop_data[dom_name] == domains[idd]),
      replace = TRUE
    )
  }

  if (sum(in_samp) != length(domains)) {
    for (idd in which(!in_samp)) {
      block_err[[idd]] <- sample(
        gb_res,
        size    = sum(pop_data[dom_name] == domains[idd]),
        replace = TRUE
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
                     cov_names = cov_names, model) {
  smp_data_tmp <- smp_data
  residuals_gb <- Y - unit_pred_smp   # ≈ e_dk (idiosyncratic residuals)

  if (is.null(model$ran_eff_sd) || model$ran_eff_sd < 1e-8) {
    # Singular LMM: all random effects are zero; use raw residuals for e*.
    gb_res_sd <- sd(residuals_gb)
    gb_res    <- if (gb_res_sd > 1e-10) (residuals_gb / gb_res_sd) * error_sd
                 else residuals_gb
    gb_res    <- gb_res - mean(gb_res)
    ran_effs  <- rep(0, length(unique(smp_data[[dom_name]])))
  } else {
    # ── ran_effs: extract actual per-domain BLUPs û_d from the fitted lme4 ──
    # ranef(effect_model)[[dom_name]] is a data.frame with one row per
    # in-sample domain and column "(Intercept)" containing û_d.
    blup_df  <- as.data.frame(lme4::ranef(model$effect_model)[[dom_name]])
    blup_vec <- setNames(blup_df[, "(Intercept)"], rownames(blup_df))

    insamp_doms <- as.character(unique(smp_data[[dom_name]]))
    ran_effs    <- as.numeric(blup_vec[insamp_doms])
    ran_effs[is.na(ran_effs)] <- 0L

    # Rescale to model$ran_eff_sd: BLUPs are shrunk toward 0 (sd(û_d) ≤ σ_u);
    # inflating them back to σ_u matches the parametric bootstrap assumption
    # that u_d* ~ distribution with sd = σ_u.
    re_sd <- sd(ran_effs)
    if (re_sd > 1e-10) ran_effs <- (ran_effs / re_sd) * model$ran_eff_sd
    ran_effs <- ran_effs - mean(ran_effs)

    # ── gb_res: idiosyncratic residuals, NOT domain-demeaned ─────────────────
    # residuals_gb = Y - unit_pred_smp = Y - gb_smp - û_d ≈ e_dk.
    # Rescale to σ_e for minor finite-sample drift.
    gb_res_sd <- sd(residuals_gb)
    gb_res    <- if (gb_res_sd > 1e-10) (residuals_gb / gb_res_sd) * error_sd
                 else residuals_gb
    gb_res    <- gb_res - mean(gb_res)
  }

  smp_data_tmp$gb_res <- gb_res
  list(gb_res = gb_res, ran_effs = ran_effs, smp_data = smp_data_tmp)
}
