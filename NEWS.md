# povmap 1.0.1

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