# The parametric-bootstrap truth generators draw one area effect per population domain
# (and one sub-area effect per population subdomain) and add it to every unit. The unit
# error must therefore be drawn at sigma_e2 for every unit. Before this fix, units in
# domains without sample observations were drawn at sigma_e2 + sigma_u2 (and, on the
# two-fold and population-weighted paths, at sigma_e2 + sigma_h2 [+ sigma_u2]) and then
# had the area draw added as well, so their variance around the regression line was
# sigma_e2 + 2 sigma_u2 and the within-domain variance in those domains was
# sigma_e2 + sigma_u2. Point estimates and in-sample MSEs were unaffected; out-of-sample
# MSEs were wrong (indicator-dependent direction).
#
# The tests reconstruct each generator's output from the RNG stream directly, so they
# check two things at once: that the fix is surgical (the number and order of draws is
# unchanged, so draws for sampled domains and the area effects are exactly what they
# were), and that the unsampled block is now scaled by sqrt(sigma_e2) alone. On the
# unfixed code the reconstruction fails on the unsampled block only.

make_pop <- function(D = 60, S = 3, n_sub = 40, seed = 1) {
  set.seed(seed)
  dom <- rep(seq_len(D), each = S * n_sub)
  sub <- rep(seq_len(D * S), each = n_sub)
  N <- length(dom)
  obs_dom <- dom <= D / 2                                  # first half of the domains sampled
  obs_subdom <- obs_dom & (((sub - 1) %% S) + 1) <= 2       # two of three subdomains observed
  list(
    N = N, D = D, S = S, dom = dom, sub = sub, mu = rnorm(N, 8, 0.5),
    obs_dom = obs_dom, obs_subdom = obs_subdom,
    fw = list(N_pop = N, obs_dom = obs_dom, N_dom_pop = D, n_pop = as.numeric(table(dom)),
              obs_subdom = obs_subdom, N_subdom_pop = D * S,
              n_pop_subdom = as.numeric(table(sub)),
              pop_domains_vec = dom, pop_subdomains_vec = sub,
              indicator_names = c("Mean", "Head_Count"), pop_weights = NULL,
              aggregate_to_vec = NULL, threshold = 8, MSE_pop_weights = "mpw",
              pop_data = data.frame(mpw = rep(4, N)))
  )
}
mp <- list(sigmae2est = 0.45, sigmau2est = 0.30, sigmah2est = 0.20)

test_that("superpopulation: unit error at sigma_e2 for every unit, draws otherwise unchanged", {
  p <- make_pop(); gm <- list(mu_fixed = p$mu)
  set.seed(123)
  r <- superpopulation(p$fw, mp, gm, lambda = NULL, shift = NULL,
                       transformation = "no", fixed = NULL)
  # the generator draws the sampled units, then the unsampled units, then the D area effects
  set.seed(123)
  z_smp <- rnorm(sum(p$obs_dom)); z_uns <- rnorm(sum(!p$obs_dom)); z_u <- rnorm(p$D)
  eps <- numeric(p$N)
  eps[p$obs_dom] <- sqrt(mp$sigmae2est) * z_smp
  eps[!p$obs_dom] <- sqrt(mp$sigmae2est) * z_uns
  vu <- sqrt(mp$sigmau2est) * z_u
  expect_equal(r$vu_tmp, vu)
  expect_equal(r$pop_income_vector, p$mu + eps + rep(vu, p$fw$n_pop))
  # stated as a variance, for the reader: within-domain variance in unsampled domains
  e <- r$pop_income_vector - p$mu - rep(vu, p$fw$n_pop)
  v_uns <- mean(tapply(e[!p$obs_dom], p$dom[!p$obs_dom], var))
  expect_lt(abs(v_uns - mp$sigmae2est), 0.05)
  expect_gt(abs(v_uns - (mp$sigmae2est + mp$sigmau2est)), 0.2)
})

test_that("superpopulation_2f: unit error at sigma_e2 in every block", {
  p <- make_pop(); gm <- list(mu_fixed = p$mu)
  set.seed(123)
  r <- superpopulation_2f(p$fw, mp, gm, lambda = NULL, shift = NULL,
                          transformation = "no", fixed = NULL)
  # draw order: D area effects, D*S subdomain effects, then the three unit blocks
  b1 <- p$obs_dom & p$obs_subdom; b2 <- p$obs_dom & !p$obs_subdom; b3 <- !p$obs_dom
  set.seed(123)
  vu <- sqrt(mp$sigmau2est) * rnorm(p$D)
  eta <- sqrt(mp$sigmah2est) * rnorm(p$D * p$S)
  eps <- numeric(p$N)
  eps[b1] <- sqrt(mp$sigmae2est) * rnorm(sum(b1))
  eps[b2] <- sqrt(mp$sigmae2est) * rnorm(sum(b2))
  eps[b3] <- sqrt(mp$sigmae2est) * rnorm(sum(b3))
  expect_equal(r$vu_tmp, vu)
  expect_equal(r$eta_tmp, eta)
  expect_equal(r$pop_income_vector,
               p$mu + eps + rep(vu, p$fw$n_pop) + rep(eta, p$fw$n_pop_subdom))
})

test_that("true_indicators_weighted: the Mean truth uses sigma_e2 for every unit", {
  p <- make_pop(); gm <- list(mu_fixed = p$mu)
  fw <- p$fw; fw$smp_subdomains <- "sub"; fw$pop_subdomains <- "sub"   # two-fold branch
  set.seed(123)
  r <- true_indicators_weighted(fw, mp, gm, lambda = NULL, shift = NULL,
                                transformation = "no", fixed = NULL)
  # draw order on this path: D*S subdomain effects, D area effects, then N unit errors with
  # a per-unit sd; rnorm(N, 0, sd) with vector sd is z * sd elementwise
  set.seed(123)
  eta <- sqrt(mp$sigmah2est) * rnorm(p$D * p$S)
  vu <- sqrt(mp$sigmau2est) * rnorm(p$D)
  eps <- sqrt(mp$sigmae2est) * rnorm(p$N)
  y <- p$mu + rep(vu, fw$n_pop) + rep(eta, fw$n_pop_subdom) + eps / sqrt(4)
  expect_equal(r$vu_tmp, vu)
  expect_equal(as.numeric(r$true_indicators[, "Mean"]),
               as.numeric(tapply(y, p$dom, mean)))
})
