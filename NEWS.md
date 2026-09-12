# povmap 1.0.1

* Parametric bootstrap MSE: the truth generator drew the unit-level error for census
  units in domains without sample observations at `sigmae2est + sigmau2est` and then
  added the domain effect, drawn at `sigmau2est`, to every unit, so those units had
  variance `sigmae2est + 2 * sigmau2est` around the regression line and their domains a
  within-domain variance of `sigmae2est + sigmau2est`. The unit error is now drawn at
  `sigmae2est` for every unit, as in Molina and Rao (2010) and in emdi's documented
  algorithm. Point estimates and in-sample MSEs are unchanged (same draws, same seed);
  out-of-sample MSEs change, in an indicator-dependent direction (on MEX2010 data the
  head-count MSE rises about 2 percent). The same construction in `superpopulation_2f()`
  and `true_indicators_weighted()` is corrected in the same way. The one-fold site is
  inherited from emdi and is in every released povmap; the other two are david3 only.
  Test: `tests/testthat/test_generator_unsampled_variance.R`.

* `xgb()`: the reported area point estimate and its bootstrap interval are now both
  computed by back-transforming each population cell and then aggregating (`hat_pc`),
  rather than aggregating on the model's transformed scale and back-transforming the
  area aggregate once (`hat`). The two orders differ whenever the transformation is
  nonlinear; for `transformation = "arcsin"` the difference is the Jensen term, which
  is positive below a rate of 0.5, negative above it, and zero at 0.5. Measured on a
  60-domain panel spanning rates 0.05-0.96 it reaches 0.069 on the rate scale and
  correlates 0.99 with the analytic leading term `cos(2*asin(sqrt(rate)))`. Because it
  changes sign at 0.5 its average over a symmetric panel is near zero, so an aggregate
  comparison will understate it. The previous value is retained as `ind$Mean_agg` (and
  `Mean_boot_agg` internally) for reconciliation. The bootstrap replicate is computed
  by reusing the already-drawn area effect, so no additional random numbers are
  consumed and the bootstrap stream is unchanged. With `transformation = "no"` the two
  orders are identical to floating point (max abs difference 2.8e-16 in test), and the
  benchmarked path already used the per-cell order, so benchmarked results are
  unaffected.

* Fit statistics: `summary.xgb`, `summary.megb`, and `xgb_cv` previously reported
  two different statistics under the single name "R2". `summary.*` reported the
  squared Pearson correlation `cor(y, yhat)^2` in-sample, while `xgb_cv` reported
  the proportion of variance explained out of sample. These are distinct and
  coincide only for in-sample OLS fits. Each function now reports BOTH statistics
  for its own setting, under unambiguous names:
    - `summary.*` `coeff_determ` gains `Squared_correlation` (the old "R2", the
      squared Pearson correlation) and `R2_prop_var` (the new `1 - SSE/SST`),
      plus the `Area_` counterparts, at both the unit and the area level.
    - `xgb_cv` keeps `r2_cv` (proportion of variance, against the noisy direct
      estimate) and gains `cor_cv` (Pearson correlation) and `cor2_cv` (its square).
  `xgb_cv`'s `r2_cv` also had a denominator inconsistency (`SSE/n` over `SST/(n-1)`),
  overstating it by a factor `(n-1)/n`; it is now exactly `1 - SSE/SST` with matched
  sum denominators. `mae_cv` is unchanged. Rank correlation and MAE were never
  ambiguous; the documentation now states their comparator (the noisy direct estimate
  out of sample). No point estimate or bootstrap interval is affected.

* `xgb()`: when `perturb_benchmark = TRUE`, the benchmark-target perturbation is now
  applied on the model's transformed (arcsin/log/poisson) scale and back-transformed,
  instead of being added on the rate scale and hard-clamped to `[0, 1]`. The rate-scale
  target SE is carried onto the transformed scale via the transform's local derivative
  (delta method), preserving the intended magnitude. This removes a clamp artefact that
  could contract bootstrap interval widths for boundary-proximate arcsin (rate)
  indicators. Non-arcsin indicators are unaffected in practice (the old log/poisson
  floor clamp never bound). Points are unchanged; only perturbed interval widths differ.

# povmap 1.0.0
  
* extension 1
* extension 2