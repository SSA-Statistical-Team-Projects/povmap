# povmap 1.0.1

* **`xgb_cv()` gains `configs`**, so cross-validation scores the configuration-averaged
  estimator: every fold's model averages the configurations as `xgb(configs = )` does
  (benchmarked once if a benchmark is passed), with the same set in every fold. With
  `configs = NULL` (default) nothing changes. Test: `tests/testthat/test_xgb_cv_configs.R`.

* **Configuration averaging for `xgb()` (opt-in).** Tuned configurations are often
  near-equivalent: their cross-validation errors differ by less than their standard
  error, yet the single best one, and with it the domain estimates, changes with the
  fold split. New argument `xgb(configs = )` takes several configurations (one per
  row, optional `weight`) and averages their domain estimates; benchmarking is applied
  once to the average. In the bootstrap each replicate uses one configuration, drawn
  with probability equal to its weight from its own stream (`config_seed`, default
  `seed`), whose fit generates the replicate's population and whose hyperparameters
  refit it, so the intervals include the choice of configuration. The replicate random
  streams are unchanged. The result adds `configs`, `config_draws`, `ind_by_config` and
  `models`. With `configs = NULL` (default) nothing changes, and a one-row `configs` is
  identical to passing the same values as arguments.
  - New `xgb_top_configs(tune, rule = "one_se_paired", floor = 3, cap = 8)` picks the
    set from one or several `xgb_tune()` results on the same grid: scored by the mean
    cross-validation error over the fold splits, kept if the mean paired difference
    from the best is within its standard error, filled to `floor` and cut to `cap`.
    `rule = "one_se"` and `"k_best"` are also available.
  - New `xgb_tune(grid = )`: score a given list of configurations (for example more
    rounds only at the smaller `eta`) instead of the full cross of the vectors.
  - In a seed check on four Nigeria ward indicators (three fold splits per replicate,
    480 configurations), the averaged estimates moved 0 to 0.4 percent of wards by more
    than half their interval width between seeds, against 13 to 19 percent for the
    single best configuration on the original 32-configuration grid.
  - Tests: `tests/testthat/test_xgb_configs.R`; the default path is covered by
    `tests/testthat/test_xgb_nga_reproduction.R`.

* **`xgb_tune()` gains `variance_y`, so tuning also mirrors `xgb()`'s
  heteroscedasticity correction.** If given, each fold model's weights are multiplied
  by `variance_y^-0.5` and divided by their mean once more, exactly as `xgb()` weights
  its fit, so the tuner and `xgb()` pass identical weights to xgboost for the same
  training rows. Pass the same `variance_y` to both. With the default `NULL` the
  weights are bit-identical to those of the `rescale_weights` fix below. Tests:
  `tests/testthat/test_xgb_tune_weights.R` (the `variance_y` case) and
  `tests/testthat/test_xgb_tune_fit_weights.R` (`min_child_weight` and `lambda`
  change the CV scores; the scores do not depend on the scale of the survey weights).

* **`xgb_tune()` now fits its candidate models on the weights `xgb()` fits on (bug fix).**
  - **The bug.** The tuner passed the raw survey weights to xgboost, while `xgb()`
    rescales them within each domain (`rescale_weights = TRUE`) and divides them by
    their mean. With weights summing to population counts, the tuner's scale was
    hundreds to thousands of times the fit's, so `min_child_weight` almost never bound
    and `lambda` barely shrank a leaf during tuning, then both acted at full strength
    in the fit. The tuned values of those two parameters were effectively arbitrary.
  - **The fix.** New argument `rescale_weights = TRUE`, mirroring `xgb()`: each fold
    model's weights are computed from its own training rows exactly as `xgb()` would
    compute them, including the division by the mean when `rescale_weights = FALSE`.
    The held-out domain means that score the candidates keep the survey weights.
    The result records `rescale_weights`.
  - **Unchanged.** The fold assignment, the random draws and `xgb()` itself, so a
    search with the same `seed` evaluates the same candidates on the same folds as
    before; only their scores change. Equal weights (all 1) give the same result as
    before.
  - **Action.** Configurations tuned with an earlier `xgb_tune()` and survey weights
    should be retuned. `tidy_xgb_tune()` is not changed.

* **`megb()` defaults change: ten domain folds, cross-fitted prediction for sampled
  domains, the fold-model average for unsampled domains, and bootstraps that predict
  their replicates the same way.**
  - **New option.** `predict_unsampled = c("foldmean", "full")`: with `"foldmean"`, an
    unsampled domain's booster part is the average of the predictions of the fold
    models of the final EM iteration, rather than the full booster's.
  - **New defaults.**
    - `cv_nfold` defaults to 10 (was 5).
    - `predict_sampled` and `predict_unsampled` default to `NULL`, which means
      `"crossfit"` and `"foldmean"` wherever fold-model prediction applies: the
      xgboost engine, `cv_folds = "domain"` and, with `mse = TRUE`, a `"leaves_only"`
      or `"full"` bootstrap. Otherwise they fall back to `"full"` with a message.
      Asking explicitly for `"crossfit"` or `"foldmean"` where they cannot apply is
      an error.
  - **The bootstraps follow the prediction.**
    - The `"leaves_only"` bootstrap refreshes every fold model and takes sampled
      domains' booster parts from their refreshed fold models and unsampled
      domains' from the average of them.
    - The `"full"` bootstrap now supports fold-model prediction: each replicate fit
      keeps its own fold models and predicts from them as the point estimate does.
  - **Why.** In a Colombian evaluation (100 replicates), the new settings lowered
    RMSE against the old defaults:
    - by 4.3 percent with 19 municipal covariates and 3.5 percent with 24 (8.4 and
      6.4 percent in sample; the fold-model average 0.8 percent out of sample);
    - neutrally with sub-area covariates.
    The fold-model bootstrap was the best calibrated in sample (bootstrap variance
    over actual MSE 1.07-1.11, against 1.25-1.27 for the five-fold full-prediction
    bootstrap).
  - **Run-time cost.** A point fit takes 1.6 to 1.9 times as long as with five folds,
    and a fit with the default bootstrap (B = 100) about 2.1 times as long as the
    five-fold, full-prediction configuration.
  - **To reproduce results from before 1.0.1,** pass `cv_nfold = 5, predict_sampled =
    "full", predict_unsampled = "full"`.
  - **What is unchanged.** The fit (random effects, variance components) does not
    depend on the prediction settings. `$crossfit$ind` reports the estimates of the
    same fit under the other settings (`Mean_full`, `Mean_crossfit`, `Mean_foldmean`).
    The internal default of `train_gbmodel()` and `em_gb_lmm()` is also ten folds.
  - **Tests:** `tests/testthat/test-megb-defaults.R`, including one confirming that,
    on unperturbed data, the bootstrap's fold-model average for unsampled domains
    reproduces the point estimate's. Earlier megb tests pin the settings they were
    written for.

* `megb()` bootstrap (`bootstrap_refit = "leaves_only"`, the default): two fixes that
  change its intervals; the point estimates are unchanged.
  - The leaf refresh now uses the survey weights (rescaled within domains) that the
    booster was trained with. It was unweighted, so even on the original data it did
    not reproduce the fitted booster.
  - Each replicate's random effects are now re-estimated the same way as the point
    estimate's: the fold models of the final EM iteration are refreshed on their own
    training rows, and the mixed model is fitted by ML to their out-of-fold residuals.
    Until now, the replicates refitted it by REML to the refreshed full booster's
    in-sample residuals. That booster has seen the replicate's own rows, so with
    domain-level covariates it can absorb the drawn domain effect.
  The fit keeps its fold models when `mse = TRUE` (about one leaf refresh per fold per
  replicate in extra cost), and works with domain or row folds. `predict_sampled =
  "crossfit"` now also has a bootstrap: each sampled domain's replicate estimate takes
  its booster part from its refreshed fold model. `"lmm_only"` and `"full"` are
  unchanged. Tests: `tests/testthat/test-megb-bootstrap-fix.R`, including one that
  confirms the bootstrap's random-effect step reproduces the point estimate's on
  unperturbed data, and one that confirms the weighted refresh reproduces the booster
  on its own data while the unweighted one does not.

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