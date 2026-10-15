# Internal helper: domain summary info used by megb_em
data_info <- function(dom_name, pop, smp) {
  in_dom    <- unique(smp[[dom_name]])
  total_dom <- unique(pop[[dom_name]])
  oosamp    <- !total_dom %in% in_dom

  list(
    n_surv  = length(smp[[dom_name]]),
    n_pop   = length(pop[[dom_name]]),
    d_out   = sum(oosamp),
    d_in    = length(in_dom),
    d_total = length(total_dom),
    ni_smp  = table(smp[[dom_name]]),
    ni_pop  = table(pop[[dom_name]]),
    dom_name = dom_name
  )
}

# Internal: conditional variance of the area random effect, per domain, from the lmer fit that
# predicts it (so the two cannot disagree). For a domain with data this is lme4's conditional
# (posterior) variance of the BLUP, sigma_b^2 (1 - gamma_i) with gamma_i = sigma_b^2 W_i /
# (sigma_b^2 W_i + sigma_e^2), W_i the domain's summed weights (n_i when unweighted); for a
# domain without data predict(allow.new.levels = TRUE) sets the effect to 0 and the variance
# is the full sigma_b^2. Returns a numeric vector aligned to `domains`.
.megb_re_condvar <- function(lmefit, domains) {
  re <- lme4::ranef(lmefit, condVar = TRUE)
  if (length(re) != 1L || ncol(re[[1L]]) != 1L)
    stop("megb: the random-effect variance correction supports a single random intercept.")
  pv  <- as.numeric(attr(re[[1L]], "postVar"))
  sb2 <- as.data.frame(lme4::VarCorr(lmefit))$vcov[1L]
  out <- rep(sb2, length(domains))
  m   <- match(as.character(domains), rownames(re[[1L]]))
  out[!is.na(m)] <- pv[m[!is.na(m)]]
  out
}

# Internal: the Jensen-corrected back-transformation used by megb() and its bootstrap. It takes the
# transformed prediction and re_var, the conditional variance of the area random effect behind it
# (see .megb_re_condvar), which is added to sigma2_e; re_var = 0 gives the unit-error-only form.
# Built with a lean captured environment (scalars + the package namespace) so parallel workers
# do not serialise the full megb() scope (~1.5 GiB).
.megb_make_corrected_bt <- function(sigma2_e, transformation, ls_lambda = NULL) {
  corrected_bt <- function(hat_t, re_var = 0) {
    s2 <- .sigma2 + re_var
    if (all(s2 == 0) || .transf %in% c("no", "logistic")) {
      return(switch(.transf,
        "no"        = hat_t,
        "log"       = exp(hat_t),
        "log.shift" = exp(hat_t) - .ls_lam,
        "arcsin"    = arcsin_transform_back(hat_t),
        "logistic"  = logit_transform_back(hat_t),
        "poisson"   = expm1(hat_t),
        hat_t
      ))
    }
    switch(.transf,
      "log"       = exp(hat_t + 0.5 * s2),
      "log.shift" = exp(hat_t + 0.5 * s2) - .ls_lam,
      "poisson"   = pmax(0, expm1(hat_t + 0.5 * s2)),
      "arcsin"    = arcsin_expected_back(hat_t, s2),
      hat_t
    )
  }
  # The random-effect variance enters only for the log, log.shift, poisson and arcsin back-transformations;
  # for 'no' and 'logistic' callers skip .megb_re_condvar (see .megb_re_var).
  attr(corrected_bt, "uses_re_var") <- transformation %in% c("log", "log.shift", "poisson", "arcsin")
  environment(corrected_bt) <- list2env(
    list(.sigma2 = sigma2_e, .transf = transformation,
         .ls_lam = if (transformation == "log.shift") ls_lambda else NULL),
    parent = getNamespace("povmap")
  )
  corrected_bt
}

# Internal: re_var argument for a corrected_bt: the conditional random-effect variances when the
# back-transformation uses them, otherwise zeros (untransformed and logistic megb run as before).
.megb_re_var <- function(corrected_bt, lmefit, domains) {
  if (isTRUE(attr(corrected_bt, "uses_re_var"))) .megb_re_condvar(lmefit, domains)
  else rep(0, length(domains))
}
