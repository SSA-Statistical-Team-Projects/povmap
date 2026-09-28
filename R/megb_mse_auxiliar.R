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

# ── Leaf-level domain aggregation for the leaves_only bootstrap ─────────────────
# (28 Sep 2026; options(povmap.leaves_only.leaf_aggregate), default TRUE.)
# The leaves_only bootstrap refreshes leaf VALUES only; every population row stays
# in the same leaf of every tree. A domain mean of the refreshed predictions is
# therefore base_score + sum over leaves of (rows of the domain in that leaf) x
# (leaf value), so it can be computed from the leaf values alone instead of
# predicting every population row in every replicate. Exact in real arithmetic;
# in floating point it differs from xgboost's single-precision accumulation at the
# 1e-7 level (see phase1_note/47 of the Colombia project for the measured sizes).

# Leaf node ids and values of every tree, plus the base score, from the JSON dump.
.megb_leaf_info <- function(booster) {
  j <- jsonlite::fromJSON(rawToChar(xgboost::xgb.save.raw(booster, raw_format = "json")),
                          simplifyVector = FALSE)
  trees <- j$learner$gradient_booster$model$trees
  list(base   = as.numeric(gsub("[][]", "", j$learner$learner_model_param$base_score)),
       leaves = lapply(trees, function(t) {
         lc <- unlist(t$left_children)
         list(leaf = which(lc == -1L) - 1L, v = unlist(t$split_conditions)[lc == -1L])
       }))
}

# Once per fit: for every (domain, leaf) pair, the number of population rows and
# the sum of their weights, plus a collapse grouping on domain that every replicate
# reuses. dom_idx (1..n_dom, every domain present) and pw must be aligned to X_pop.
.megb_leaf_aggregator <- function(booster, X_pop, dom_idx, n_dom, pw = NULL, chunk = 40L) {
  info  <- .megb_leaf_info(booster)
  nt    <- length(info$leaves)
  nleaf <- vapply(info$leaves, function(z) length(z$leaf), 1L)
  off   <- c(0L, cumsum(nleaf))
  dm    <- xgboost::xgb.DMatrix(X_pop)
  w     <- if (is.null(pw)) rep(1, nrow(X_pop)) else as.numeric(pw)
  pd <- pg <- pn <- pw_ <- list()
  for (a in seq(1L, nt, by = chunk)) {
    b  <- min(nt, a + chunk - 1L)
    lf <- predict(xgboost::xgb.slice.Booster(booster, a, b), dm, predleaf = TRUE)
    if (is.null(dim(lf))) lf <- matrix(lf, ncol = 1L)
    g  <- vapply(seq_len(ncol(lf)), function(k) {
      t <- a + k - 1L; match(lf[, k], info$leaves[[t]]$leaf) + off[t] }, numeric(nrow(lf)))
    if (anyNA(g)) stop("leaf aggregation: a population row fell in a node that is not a leaf")
    G  <- collapse::GRP(list(d = rep(dom_idx, ncol(lf)), g = as.integer(g)), sort = FALSE)
    k  <- length(pd) + 1L
    pd[[k]]  <- G$groups$d
    pg[[k]]  <- G$groups$g
    pn[[k]]  <- collapse::GRPN(G, expand = FALSE)
    pw_[[k]] <- collapse::fsum(rep(w, ncol(lf)), G, use.g.names = FALSE)
  }
  d <- unlist(pd, use.names = FALSE)
  list(g = unlist(pg, use.names = FALSE), n = unlist(pn, use.names = FALSE),
       w = unlist(pw_, use.names = FALSE), GD = collapse::GRP(d),   # sorted: groups 1..n_dom
       n_dom = n_dom, n_leaves = off[nt + 1L],
       nd = tabulate(dom_idx, n_dom), sw = collapse::fsum(w, dom_idx, use.g.names = FALSE),
       weighted = !is.null(pw))
}

# Per replicate: unweighted and weighted domain means of a refreshed booster.
.megb_leaf_domain_means <- function(LA, booster) {
  info <- .megb_leaf_info(booster)
  v <- unlist(lapply(info$leaves, `[[`, "v"))
  if (length(v) != LA$n_leaves) stop("leaf aggregation: the refreshed booster's tree structure changed")
  vg  <- v[LA$g]
  s_n <- collapse::fsum(vg * LA$n, LA$GD, use.g.names = FALSE)
  s_w <- collapse::fsum(vg * LA$w, LA$GD, use.g.names = FALSE)
  list(unweighted = info$base + s_n / LA$nd, weighted = info$base + s_w / LA$sw)
}
