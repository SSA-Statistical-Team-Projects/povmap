
#' MSE helper functions
#' @param domains domain
#' @param in_samp in_samp
#' @param smp_data sample data
#' @param dom_name domain name
#' @param pop_data population data
#' @param gb_res gb_res
#' @importFrom dplyr left_join group_by summarise
#' @importFrom rlang sym


# Block-sampling errors in the MSE procedures ----------------------------------------------
block_sample <- function(domains, in_samp, smp_data, dom_name, pop_data, gb_res) {
  block_err <- vector(mode = "list", length = length(domains))

  for (idd in which(in_samp)) {
    block_err[[idd]] <- sample(
      gb_res[smp_data[dom_name] == domains[idd]],
      size = sum(pop_data[dom_name] == domains[idd]),
      replace = TRUE
    )
  }

  if (sum(in_samp) != length(domains)) {
    for (idd in which(!in_samp)) {
      block_err[[idd]] <- sample(
        gb_res,
        size = sum(pop_data[dom_name] == domains[idd]),
        replace = TRUE
      )
    }
  }

  unlist(block_err)
}

# Draws random survey samples in the MSE procedures ----------------------------------------
sample_select <- function(pop, smp, dom_name) {
  smpSizes <- table(smp[dom_name])
  smpSizes <- data.frame(
    smpidD = as.character(names(smpSizes)),
    n_smp = as.numeric(smpSizes),
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
    do.call(
      rbind,
      mapply(dfList, ns, FUN = function(df, n) {
        popInd <- seq_len(nrow(df))
        sel <- base::sample(popInd, n, replace = FALSE)
        df[sel, ]
      }, SIMPLIFY = FALSE)
    )
  }

  samples <- stratSamp(splitPop, smpSizes$n_smp)
  samples
}

# Computes REB random components in the MSE procedures ------------------------------------
ran_comp <- function(unit_pred_smp, smp_data, Y, dom_name, error_sd, cov_names = cov_names, model) {
  smp_data_tmp <- smp_data
  residuals_gb <- Y - unit_pred_smp
  smp_data_tmp$gb_res <- residuals_gb

  # Random Effects unscaled
  ran_effs_unscaled <- smp_data_tmp %>%
    group_by(!!sym(dom_name)) %>% # Group by domain variable
    dplyr::summarise(r_bar = mean(gb_res, na.rm = TRUE)) # Compute mean

  smp_data_res <- dplyr::left_join(smp_data_tmp, ran_effs_unscaled, by = dom_name)
  smp_data_res$gb_eij <- smp_data_res$gb_res - smp_data_res$r_bar

  # Standardize gb_res
  gb_res <- smp_data_res$gb_eij
  gb_res <- (gb_res / sd(gb_res)) * error_sd
  gb_res <- gb_res - mean(gb_res) # Center

  # Standardize random effects
  ran_effs <- ran_effs_unscaled$r_bar
  ran_effs <- (ran_effs / sd(ran_effs)) * model$ran_eff_sd
  ran_effs <- ran_effs - mean(ran_effs) # Center

  list(
    gb_res = gb_res,
    ran_effs = ran_effs,
    smp_data = smp_data_res
  )
}

