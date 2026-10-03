# povmap 1.0.1

* `megb()`: new argument `cv_nfold` (default 5), the number of folds of the internal
  cross-validation (xgboost), of whole domains or of rows according to `cv_folds`. It is
  passed through to the fit and to the `bootstrap_refit = "full"` refits. With the default,
  `xgb.cv` receives exactly the historical arguments and fold assignment. More folds train
  each fold model on a larger share of the sample, at roughly proportional cost per EM
  iteration. Under `predict_sampled = "crossfit"` there is one fold model per fold, and ten
  folds removed the small in-sample cost cross-fitting had with five in a Colombian
  evaluation, while keeping its gain with area-level covariates. When omitted, the default
  comes from `getOption("megb.cv_nfold", 5L)`. Test: `tests/testthat/test-megb-crossfit.R`.

* `megb()`: new argument `predict_sampled = c("full", "crossfit")`. `"full"` (default)
  is unchanged. Under `"crossfit"`, each sampled domain's population rows are predicted by
  the fold model of the final EM iteration that held that domain out, plus its random
  effect. Those are the fold models whose out-of-fold residuals the random effects were
  estimated from, so the booster part and the random effect share one baseline. Under
  `"full"`, the booster has seen the domain's own rows; with covariates constant within
  domains (area-level covariates) it can learn the domain's level from them, and part of
  the deviation is then counted twice. Unsampled domains keep the full booster.
  `$crossfit$ind` also returns `Mean_full` (the `"full"` estimates from the same fit) and
  `Mean_foldmean` (unsampled domains predicted by the mean of the five fold models). The
  fit itself, including the random effects and variance components, is the same under
  both settings. **`"crossfit"` requires `cv_folds = "domain"`** (with row folds, no fold
  model excludes a whole domain), the xgboost engine, and `mse = FALSE`, since there is no
  bootstrap consistent with it yet; any other combination is an error. In a Colombian
  evaluation (20 replicates), it reduced in-sample RMSE by 7.3 percent with 19
  municipal-level covariates. With sub-area covariates it made no difference (+0.4
  percent, within about one or two standard errors). Test:
  `tests/testthat/test-megb-crossfit.R`.

* `estimators()`, `write.excel()` and `write.ods()` with `CV = TRUE` (or `MSE`/`var = TRUE`)
  now pair each variance or MSE column with the point estimate it measures by name, not by
  column position. Since the `xgb()` change below, `$ind` carries `Mean_agg` while `$var`
  has no variance for it, so the CV division stopped with "'/' only defined for
  equally-sized data frames". Even had the division gone through, the positional pairing
  would have labelled `Var_bench` as `Mean_agg_Var`. `$var$Mean` is paired with `$ind$Mean`,
  and `$var$Var_bench` with `$ind$Mean_bench`, so `Mean_bench_CV` is
  `sqrt(Var_bench) / Mean_bench`, as before. `Mean_agg` is exported without `_Var` or `_CV`
  columns. Rows are matched on `Domain`. A variance column with no matching point estimate
  is now an error. Objects without `Mean_agg` export exactly as before. Estimation is
  unchanged. Test: `tests/testthat/test_cv_export_pairing.R`.

* `megb()`: new arguments `cv_folds` and `early_stopping` control the internal
  cross-validation (`xgboost::xgb.cv`) whose out-of-fold residuals feed the random-effect
  and error variances. **The default changes to `cv_folds = "domain"`**, which holds out
  5 folds of whole domains; `"rows"` restores the previous 5 folds of random rows.
  `early_stopping = TRUE` (default) is unchanged. With row folds and covariates that are
  constant within domains, each held-out row's domain is in the training folds, the
  out-of-fold residuals lose their between-domain variation, and the random-effect
  variance is estimated at or near zero, which makes the bootstrap intervals far too
  narrow (in a Colombian evaluation, zero in about half the fits and out-of-sample
  coverage of 0.44-0.53 for 95 percent intervals; about 0.88 with domain folds, at the
  same point accuracy). When the arguments are omitted their defaults come from
  `getOption("megb.cv_folds", "domain")` and `getOption("megb.early_stopping", TRUE)`,
  so code that set those options still works. If domain folds are requested but no
  domain labels reach the booster, row folds are used with a message. The settings are
  passed explicitly to the `bootstrap_refit = "full"` refits, rather than relying on
  global options reaching parallel workers. Tests: `tests/testthat/test-megb-cv-options.R`.

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