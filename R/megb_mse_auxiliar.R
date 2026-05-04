# Internal MSE bootstrap helpers for megb_em
#' @importFrom dplyr left_join group_by summarise

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
ran_comp <- function(unit_pred_smp, smp_data, Y, dom_name, error_sd,
                     cov_names = cov_names, model) {
  smp_data_tmp      <- smp_data
  residuals_gb      <- Y - unit_pred_smp
  smp_data_tmp$gb_res <- residuals_gb

  ran_effs_unscaled <- smp_data_tmp |>
    dplyr::group_by(.data[[dom_name]]) |>
    dplyr::summarise(r_bar = mean(gb_res, na.rm = TRUE))

  smp_data_res       <- dplyr::left_join(smp_data_tmp, ran_effs_unscaled, by = dom_name)
  smp_data_res$gb_eij <- smp_data_res$gb_res - smp_data_res$r_bar

  if (is.null(model$ran_eff_sd) || model$ran_eff_sd < 1e-8) {
    # Singular LMM fit: skip demeaning so block-sampling preserves domain-level
    # variance that would otherwise be lost when ran_effs are all zero.
    gb_res   <- residuals_gb
    gb_res   <- (gb_res / sd(gb_res)) * error_sd
    gb_res   <- gb_res - mean(gb_res)
    ran_effs <- rep(0, nrow(ran_effs_unscaled))
  } else {
    gb_res  <- smp_data_res$gb_eij
    gb_res  <- (gb_res / sd(gb_res)) * error_sd
    gb_res  <- gb_res - mean(gb_res)

    ran_effs <- ran_effs_unscaled$r_bar
    ran_effs <- (ran_effs / sd(ran_effs)) * model$ran_eff_sd
    ran_effs <- ran_effs - mean(ran_effs)
  }

  list(gb_res = gb_res, ran_effs = ran_effs, smp_data = smp_data_res)
}
