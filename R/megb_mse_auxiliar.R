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


# The fit's fold models for the leaves_only bootstrap (NULL when it kept none).
.megb_fold_fit <- function(model) {
  if (is.null(model$fold_models)) return(NULL)
  list(fold_models = model$fold_models, fold_rows = model$fold_rows,
       fold_domains = model$fold_domains, oof_prediction = model$oof_prediction)
}

# Leaf refresh of an xgboost booster on label y (tree structure unchanged), with
# observation weights w (the rescaled survey weights the booster was trained with)
# or unweighted when w is NULL. With the booster's own training data, label and
# weights it returns the booster unchanged.
.megb_refresh <- function(bst, X, y, w, nr) {
  dt <- xgboost::xgb.DMatrix(X, label = y)
  if (!is.null(w)) xgboost::setinfo(dt, "weight", as.numeric(w))
  suppressMessages(suppressWarnings(xgboost::xgb.train(
    params = list(updater = "refresh", process_type = "update", refresh_leaf = 1,
                  objective = "reg:squarederror"),
    data = dt, nrounds = nr, xgb_model = bst, verbose = 0)))
}

# The point estimate's random-effect step (em_gb_lmm): ML, survey-weighted, on
# out-of-fold residuals r.
.megb_oof_re_fit <- function(fit_data, lmm_formula, r, w) {
  fit_data$r <- r
  .w <- if (!is.null(w)) as.numeric(w) else NULL
  # lme4 evaluates `weights` in the data, then in the formula's environment
  environment(lmm_formula) <- environment()
  suppressMessages(suppressWarnings(lme4::lmer(
    lmm_formula, data = fit_data, REML = FALSE, weights = .w)))
}

# leaves_only bootstrap setup, once per fit. sort_smp maps the bootstrap's
# domain-sorted sample rows to the fit's row order (sorted row p is fit row
# sort_smp[p]). Each fold model is sliced to the boosting rounds its predictions
# use (xgb.cv trains fold models past the best round under early stopping, and
# predict() stops there) and must reproduce the fit's out-of-fold predictions of
# its held-out rows exactly. For crossfit prediction, also the population rows and
# positions (in `domains`) of its domains and, for the identity back-transform, a
# leaf aggregator over those rows.
.megb_fold_boot_setup <- function(fold_fit, X_smp, sort_smp, d_smp, X_pop, pop_dom_idx,
                                  domains, pw_pop, use_leaf, predict = FALSE,
                                  predict_unsampled = FALSE) {
  fm <- fold_fit$fold_models; fr <- fold_fit$fold_rows
  if (is.null(fm) || length(fm) != length(fr)) stop("leaves_only bootstrap: no fold models.")
  n <- nrow(X_smp)
  inv <- match(seq_len(n), sort_smp)                     # fit row j -> sorted row inv[j]
  oof_sorted <- fold_fit$oof_prediction[sort_smp]
  if (length(oof_sorted) != n || anyNA(inv)) stop("leaves_only bootstrap: fold rows do not match the sample.")
  if (predict && is.null(fold_fit$fold_domains))
    stop("cross-fitted bootstrap needs folds of whole domains.")
  ## unsampled domains (positions in `domains`) and their population rows, for the
  ## fold-model average (predict_unsampled)
  posU <- if (predict_unsampled && !is.null(pop_dom_idx)) which(!(domains %in% unique(d_smp))) else integer(0)
  ipU  <- if (length(posU)) which(pop_dom_idx %in% posU) else integer(0)
  locU <- if (length(posU)) match(pop_dom_idx[ipU], posU) else integer(0)
  folds <- lapply(seq_along(fm), function(k) {
    b  <- fm[[k]]
    bi <- xgboost::xgb.attr(b, "best_iteration")
    nr <- if (is.null(bi)) xgboost::xgb.get.num.boosted.rounds(b) else as.integer(bi) + 1L
    sl <- xgboost::xgb.slice.Booster(b, 1L, nr)
    ho <- sort(inv[fr[[k]]])
    if (!identical(as.numeric(predict(sl, X_smp[ho, , drop = FALSE])), as.numeric(oof_sorted[ho])))
      stop("leaves_only bootstrap: fold model ", k, " does not reproduce the fit's out-of-fold predictions.")
    out <- list(bst = sl, nr = nr, tr = setdiff(seq_len(n), ho), ho = ho)
    if (predict) {
      fd  <- fold_fit$fold_domains[[k]]
      if (!setequal(unique(d_smp[ho]), fd)) stop("cross-fitted bootstrap: fold ", k, " rows and domains disagree.")
      out$pos <- match(fd, domains); out$ip <- which(pop_dom_idx %in% out$pos)
      if (use_leaf) {
        loc <- match(pop_dom_idx[out$ip], out$pos)
        out$LA <- .megb_leaf_aggregator(sl, X_pop[out$ip, , drop = FALSE], loc, length(out$pos),
                                        pw = if (is.null(pw_pop)) NULL else pw_pop[out$ip])
      }
    }
    if (length(posU) && use_leaf)
      out$LAu <- .megb_leaf_aggregator(sl, X_pop[ipU, , drop = FALSE], locU, length(posU),
                                       pw = if (is.null(pw_pop)) NULL else pw_pop[ipU])
    out
  })
  if (anyDuplicated(unlist(lapply(folds, `[[`, "ho"))) || length(unlist(lapply(folds, `[[`, "ho"))) != n)
    stop("leaves_only bootstrap: the folds do not partition the sample.")
  list(folds = folds, pos_all = unlist(lapply(folds, `[[`, "pos")),
       ip_all = unlist(lapply(folds, `[[`, "ip")), posU = posU, ipU = ipU)
}

# Per replicate: refresh every fold model on its training rows with weights w
# (refresh = FALSE uses the fold models as fitted), predict its held-out rows
# (out-of-fold) and, with predict, its domains' population means (leaf
# aggregation) or rows, in the order of CFB$pos_all / CFB$ip_all.
.megb_fold_boot_step <- function(CFB, X_smp, X_pop, y, w, refresh = TRUE, use_leaf = TRUE,
                                 predict = FALSE, predict_unsampled = FALSE) {
  oof <- rep(NA_real_, length(y)); un <- wt <- gp <- list()
  .u <- predict_unsampled && length(CFB$posU)
  uu <- uw <- if (.u) numeric(length(CFB$posU)) else NULL
  ug <- if (.u && !use_leaf) numeric(length(CFB$ipU)) else NULL
  for (k in seq_along(CFB$folds)) {
    f  <- CFB$folds[[k]]
    fb <- if (refresh) .megb_refresh(f$bst, X_smp[f$tr, , drop = FALSE], y[f$tr],
                                     if (is.null(w)) NULL else w[f$tr], f$nr) else f$bst
    oof[f$ho] <- predict(fb, X_smp[f$ho, , drop = FALSE])
    if (.u) {                                    # unsampled: accumulate for the fold-model average
      if (use_leaf) { LU <- .megb_leaf_domain_means(f$LAu, fb); uu <- uu + LU$unweighted; uw <- uw + LU$weighted }
      else ug <- ug + as.numeric(predict(fb, X_pop[CFB$ipU, , drop = FALSE]))
    }
    if (!predict) next
    if (use_leaf) {
      LM <- .megb_leaf_domain_means(f$LA, fb); un[[k]] <- LM$unweighted; wt[[k]] <- LM$weighted
    } else gp[[k]] <- as.numeric(predict(fb, X_pop[f$ip, , drop = FALSE]))
  }
  if (anyNA(oof)) stop("leaves_only bootstrap: a sample row was held out by no fold model.")
  K <- length(CFB$folds)
  list(oof = oof, unweighted = unlist(un), weighted = unlist(wt), gb_pop = unlist(gp),
       u_unweighted = if (.u && use_leaf) uu / K else NULL, u_weighted = if (.u && use_leaf) uw / K else NULL,
       u_gb_pop = if (.u && !use_leaf) ug / K else NULL)
}

# Interval from the bootstrap prediction errors (5 Oct 2026). `err` is D x B:
# each replicate's estimate minus its bootstrap truth, e* = tau_b - tau*. The
# errors are centred per domain (the bootstrap bias is not carried into the
# interval) and their variance is returned as the domain's variance. Treating e*
# as the law of the error e = theta_hat - theta, the interval for theta is
#   [theta_hat - q(1 - alpha/2), theta_hat - q(alpha/2)],
# the basic bootstrap interval: the error quantiles are SUBTRACTED. Until
# 5 Oct 2026 megb added them ([theta_hat + q(alpha/2), theta_hat + q(1 - alpha/2)]),
# which is the same interval only when the errors are symmetric; with skewed
# errors it put the long tail on the wrong side.
.megb_error_interval <- function(err, alpha) {
  err <- err - rowMeans(err, na.rm = TRUE)
  list(var  = apply(err, 1, stats::var, na.rm = TRUE),
       q_lo = apply(err, 1, stats::quantile, probs = alpha / 2,     na.rm = TRUE, names = FALSE),
       q_hi = apply(err, 1, stats::quantile, probs = 1 - alpha / 2, na.rm = TRUE, names = FALSE))
}

# The leaf refresh for one leaves_only replicate, with every tree fixed
# (5 Oct 2026). The default (max_iterations = 0) is the single pass; with
# options(povmap.leaves_only.em = TRUE) it is the point fit's EM. y: the replicate outcome y*, in the sorted sample order. Mirrors
# em_gb_lmm: a random-effect-only ML fit on y* gives the first target, y* minus
# the estimated random effect; then each iteration refreshes the fold models on
# the target (.megb_fold_boot_step, each on its own training rows), fits the
# random effect by ML on the out-of-fold residuals y* - oof, and forms the next
# target y* - (fitted - fixed intercept). It stops when the relative change in
# the log-likelihood is at most error_tolerance or after max_iterations, as
# em_gb_lmm does. max_iterations = 0: the single pass, the default (fold models
# refreshed on y* itself, one random-effect fit). Returns the last iteration's
# fold step (cfr), random-effect fit (lmer), the target its fold models were
# refreshed on, and the number of iterations.
.megb_leaf_em <- function(CFB, X_smp, X_pop, y, w_ref, fit_data, lmm_formula, w_lmm,
                          use_leaf, predict, predict_unsampled, max_iterations, error_tolerance,
                          start_re = NULL) {
  # start_re (opt-in, options(povmap.leaves_only.em_start = "point"), EM only): start the EM from
  # the point fit's random effect per sample row instead of a random-effect-only fit on y*.
  step <- function(target) .megb_fold_boot_step(CFB, X_smp, X_pop, target, w_ref, refresh = TRUE,
                                                 use_leaf = use_leaf, predict = predict,
                                                 predict_unsampled = predict_unsampled)
  if (max_iterations < 1L) {
    cfr <- step(y)
    return(list(cfr = cfr, lmer = .megb_oof_re_fit(fit_data, lmm_formula, y - cfr$oof, w_lmm),
                target = y, iterations = 0L))
  }
  re_part <- function(m) as.numeric(stats::fitted(m)) - as.numeric(lme4::fixef(m)[1])
  target <- if (!is.null(start_re)) y - start_re else y - re_part(.megb_oof_re_fit(fit_data, lmm_formula, y, w_lmm))
  old_ll <- 0; it <- 0L
  repeat {
    it  <- it + 1L
    cfr <- step(target)
    lm1 <- .megb_oof_re_fit(fit_data, lmm_formula, y - cfr$oof, w_lmm)
    ll  <- as.numeric(stats::logLik(lm1))
    go  <- abs((ll - old_ll) / old_ll) > error_tolerance && it < max_iterations
    old_ll <- ll
    if (!go) break
    target <- y - re_part(lm1)
  }
  list(cfr = cfr, lmer = lm1, target = target, iterations = it)
}
