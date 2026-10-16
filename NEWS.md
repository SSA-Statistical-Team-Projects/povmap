# povmap 1.0.4 (v0.6.2)

* **`megb()` can perturb the benchmark target (`perturb_benchmark`, `benchmark_target_se`, `allow_zero_benchmark_se`).** The same method as `xgb()`: after the bootstrap, a draw from N(0, SE_g^2) for each benchmark group is added to the benchmarked replicate of every domain in the group, on the model's transformed scale (increment = draw times `.perturb_slope()`, back-transformed; added directly for `transformation = "no"`). The truth side is not perturbed, so the benchmarked intervals carry the sampling error of the targets. The B x G draws are made once from `seed` on a stream of their own (`seed + 7919`) and the caller's RNG state is restored, so the bootstrap itself is the same with and without perturbation. `benchmark_target_se` (a named vector keyed by `benchmark_level` group, validated as in `xgb()`) is required; there is no internal Horvitz-Thompson fallback. The result then also holds `Var_bench_unperturbed` and `benchmark_perturbation`. With the default `perturb_benchmark = FALSE` the output is identical to 1.0.3 (new fixture `megb_bench_identity_v061.rds`: benchmarked `ind`, `var`, `CI`, `mse_type` "var" and "mse", arcsin/logit_raking, log/ratio, no/ratio).
* **`megb()` takes `em_maxit` (default 25) and `em_tol` (default 1e-04).** The EM iteration cap and the stopping tolerance on the relative log-likelihood change were package constants (`.megb_em_max_iterations`, `.megb_em_error_tolerance`); they are now arguments, passed to the point fit, to the `bootstrap_refit = "full"` replicate fits and to the bootstrap's EM settings. The defaults are the old constants, so the output is unchanged.

# povmap 1.0.3 (v0.6.1 fixes on top of the convention release)

* **Arcsin back-transformation (`ebp`, `megb`).** The mean of `sin(x)^2` under normal error is now the exact `1/2 - 1/2 cos(2 mu) exp(-2 v)` (new internal `arcsin_expected_back()`), replacing a 9th-order series that had the wrong sign on the `v^2` term and lacked the double-factorial factors on the higher terms. Arcsin estimates move; the largest effect is at large `v`.
* **`ebp()` at `L > 0` runs again.** The `Quantiles` indicator returned a different number of values from the indicator-name list since d5c0068/b014b53 (9 probabilities, 8 names, and a mis-indexed `keepthese`), so `monte_carlo()` failed on a replacement-length error. `Quantiles` now returns the 7 values `Quantile_5` ... `Quantile_95` (Minimum and Maximum are their own indicators), the stray `"Quantile"` name is gone, and `monte_carlo()` sizes its results array from what the indicator functions return and stops with a clear message if that disagrees with the names.
* `sae` is in Suggests (the `hdp` example needs it).
* **`megb()` back-transformation adds the random effect's conditional variance (log, log.shift, poisson, arcsin).** The Jensen correction used only the unit-error variance `sigma_e^2`. It now uses `sigma_e^2 + v_d`, where `v_d` is the conditional variance of the area effect from the same `lmer` fit that predicts it (`lme4::ranef(condVar = TRUE)`, so weights are included): `sigma_b^2` for unsampled domains (effect set to 0) and `sigma_b^2 (1 - gamma_d)` for sampled ones. The bootstrap MSE applies it to each replicate's estimate from the replicate's own fit. Estimates for unsampled domains rise by `exp(sigma_b^2 / 2)` on the log scale. The arcsin mean also clamps `mu` to `[0, pi/2]`, as the plain back-transform does.
* **`megb()` bootstrap draws `y_star` take the plain back-transform.** The drawn `y_star` already contains the unit error `e*`, like the observed `y`, but the benchmark target built from it was also Jensen-corrected, which counted `sigma_e^2` twice in the bootstrap samples of transformed `megb` and affected its MSE. The target now back-transforms `y_star` without the correction. The bootstrap truth is unchanged: it back-transforms the drawn area effect with `sigma_e^2` only. The random-effect variance is computed only for the log, log.shift, poisson and arcsin transformations; untransformed (`"no"`) and logistic `megb` run as before (their `ind` and MSE are identical to 1.0.2).
* **Known gap: `megb(transformation = "logistic")`.** It back-transforms with the plain `plogis` and no variance correction (neither `sigma_e^2` nor the random effect's conditional variance), so estimates are biased away from 0.5: too low where p < 0.5, too high where p > 0.5. Its MSE uses a delta-method variance on the transformed-scale mean. Not changed in this release. If fixed later, use Gauss-Hermite quadrature over N(mu, sigma_e^2 + v_d), matching the exact arcsin treatment, not a probit approximation.

# povmap 1.0.2 (shared XGBoost estimation convention v0.6)

* **Defaults now follow the convention** (the arguments are unchanged; only their defaults): `xgb(area_pool = "out_of_fold", weightedBS = FALSE)`
  (`sub_pool = "in_sample"` and `area_deflate = FALSE` were already the defaults); `xgb_top_configs(floor = 1)`, so the paired one-SE set is not padded
  (`cap = 8`); `xgb_tune(n_iter = 150, nrounds_max = 3000)` and the default random-search bounds eta 0.01-0.3, max_depth 1-12, subsample 0.5-1,
  colsample_bytree 0.1-1, min_child_weight 1-1000, lambda 0.01-10000, nround 50-3000. `xgb_tune(search = )` stays `"grid"`, so existing grid calls are unchanged.
* **New `xgb_tune_repeated()`**: one draw of `n_iter` configurations (`draw_seed`), scored on fold seeds 123, 124, 125; the result goes to `xgb_top_configs()`.
* Calls that relied on the old defaults (`area_pool`, `weightedBS`, `floor`, the random-search bounds) give different results; pass the old value to reproduce them.

# povmap 1.0.1

* **`xgb()`: out-of-fold residual pools (opt-in).** The bootstrap population draws from two residual pools, the area pool
  (domain-mean residuals) and the sub-area pool (unit residual less domain mean). They were always the fit's own
  residuals; with `area_pool = "out_of_fold"` and/or `sub_pool = "out_of_fold"` they come from models that did not train on
  the unit. Area folds are grouped by domain (`pool_folds`, default 10, seeded by `pool_seed`, balanced; a data.frame
  gives an assignment); sub-area folds are grouped by sub-area within domains (a domain's sub-areas get different folds,
  so a held-out sub-area's model was trained on its siblings). `sub_pool = "out_of_fold_domain"` takes the sub-area pool
  from the domain-grouped fits (the DRC scratch-hook construction; no extra fits). `area_deflate`/`sub_deflate`
  remove the sampling error: pool values are multiplied by `sqrt(max(msr - msv, 0) / msr)`, msr the pool's mean
  squared residual and msv its mean sampling variance (bootstrap draw weights, model scale). **Sampling variances**
  (`area_sampling_var`, `sub_sampling_var`) are on the RESPONSE scale (squared standard errors) and are moved to the model
  scale by the delta method with the package slope `.perturb_slope`. `area_sampling_var`: a numeric vector NAMED by the
  domain labels (as in `smp_data[[domains]]`) holding the variance of each domain's direct (weighted-mean) estimate, or the name
  of a `smp_data` column constant within domain; domains of the pool without an entry are left out of both means (warning);
  `NULL` computes the cluster-robust (Taylor, with-replacement) variance of the weighted mean with the sample rows as the clusters.
  `sub_sampling_var`: the sampling variance of each sub-area MEAN, a numeric vector with one element per row of `smp_data` or a column name, REQUIRED for `sub_deflate` (`variance_y` is a relative heteroscedasticity weight, 1/n in DRC, not a sampling variance, and is not used). Domains without a variance (one sampled EA, or absent from the vector) are imputed by default (`impute_single_ea = TRUE`): the pooled within-domain EA variance / the domain's EA count; the number imputed is returned and reported; `FALSE` drops them with a warning. The deflation means use the weights the bootstrap draws the pool with.
  With `configs`, each configuration has its own pools and a replicate uses those of the configuration it draws. The pools are
  built after the point estimate and with the random-number state saved and restored, so **point estimates do not depend
  on the pool and an existing call is unchanged** (checked: `identical()` on the whole estimates object, pre-change build
  vs this one, for log, arcsin and poisson, with and without benchmarking, perturbation, `variance_y` and `configs`). The folds,
  pools and factors are returned in `$pools` (only when an option is used).
  - `area_pool = "by_status"`: in the bootstrap every census domain draws one area residual from the pool (a sampled domain's own residual is not attached
    to it); under `by_status` a sampled domain draws from the in-sample pool and a census-only domain from the out-of-fold pool, at the same pool position
    (same domains, same draw probabilities), so the random-number stream is that of the other settings. Without a benchmark the sampled domains' results equal
    `in_sample`'s and the census-only domains' equal `out_of_fold`'s (tested at a fixed seed); with a benchmark the in-bootstrap re-estimation and the raking
    mix the domains of a province, so only the draws are guaranteed equal. With `area_deflate`, each pool is deflated separately with its own residuals,
    mean squared residual and the bootstrap's draw weights: `$pools$configs[[k]]$area_deflation` (in-sample pool) and `area_deflation_alt` (out-of-fold pool),
    each with factor, msr, msv, n_missing and n_imputed. Works with `configs` and `pool_tune`.
  - `deflate_from = "residuals"` (default `"outcomes"`, unchanged): the sampling variance of the area-pool deflation is the cluster-robust variance of the
    domain's weighted mean of the unit residuals of the fit that produced the pool, on the model (transformed) scale, no conversion:
    `m/(m-1) * sum_c (sum_{i in c} w_i (r_i - rbar))^2 / (sum w)^2`, clusters = sample rows grouped by `sub_domains`. A one-cluster domain is imputed (`impute_single_ea`:
    pooled within-domain unit-residual variance per cluster / the domain's cluster count; the count is returned). Area pool only: the sub-area route
    (`sub_sampling_var`) is unchanged. Not allowed together with `area_sampling_var`.
  - `xgb_tune(keep_oof = TRUE)` keeps each candidate's out-of-fold predictions; `xgb(pool_tune = tune, area_pool = "out_of_fold")`
    builds the area pool from them with no new fits (saving `pool_folds` fits per configuration), after checking that the tune
    is the same model on the same data (transformation, weights, `variance_y`, covariates, domains, outcome, hyperparameters and
    rounds, and sample covariates equal to the census cells') and stopping with the reason otherwise.
  - Test: `tests/testthat/test_xgb_oof_pools.R`; the DRC scratch-hook construction (`oof_lib.R`) is reproduced exactly.

* **Least-overfit configuration selection (opt-in).** `xgb_top_configs()` gains `rule = "least_overfit"`:
  among the configurations of the paired one-SE set (floor 1, no cap), the single one with the lowest
  out-of-fold to in-sample residual ratio, ties to the better score; and `rule = "one_se_ratio"`: the
  set members with ratio at most `max_ratio` (default 2), averaged with equal weights (at most `cap`),
  falling back to the single lowest-ratio configuration when fewer than two qualify. New arguments
  `ratio`, `max_ratio` and `ratio_cap` (the best 40 set members are ratioed, as in the DRC selection);
  the ratio rules add an `overfit_ratio` column that `xgb(configs = )` and `xgb_cv(configs = )` ignore.
  The ratio comes from the new exported `xgb_overfit_ratio()` (full fit vs territory-grouped
  cross-validation, survey-weighted residual SD on the transformed scale, total, domain and
  within-domain pools) or from `xgb_tune(overfit_ratio = TRUE)`, which adds an `overfit` element on the
  search's own folds. **Nothing changes for existing callers**: `overfit_ratio` defaults to `FALSE`, the
  default `rule` is still `"one_se_paired"`, and outputs of the existing rules are identical to before.
  - **Residual pools.** `xgb_overfit_ratio()` returns all three pools in one object (printed together):
    `total` (unit-level residuals, the pool the DRC selection ranks on), `domain` (domain-mean
    residuals) and `within_domain` (unit residual less domain mean). `xgb_top_configs(ratio_pool = )`
    says which one is ranked and thresholded (default `"total"`); quote the pool with the number.
  - `max_ratio = 2` is a provisional, unvalidated default. If out-of-fold residual pools are adopted the
    threshold only needs to exclude badly overfitted configurations; observed ratios on comparable
    pools run about 1.1 to 1.9 (DRC, total pool, lowest-ratio member of each set).
  - `xgb_tune()` gains two opt-in arguments, both inert by default: `score_scale = "response"` scores the
    held-out domain means on the outcome's own scale (back-transformed predictions; with it,
    `transformation = "poisson"` trains on `log1p(y)` as `xgb()` does, whereas `xgb_tune()` otherwise
    leaves a poisson outcome untransformed), and `fit_seed` fixes the xgboost seed separately from the
    fold seed `seed`.
  - `inst/scripts/least_overfit_rescore.R` re-runs the DRC selection from the saved candidate tables with
    package code only (fresh folds 901 and 902, ratios, `rule = "least_overfit"`).
  - Test: `tests/testthat/test_xgb_least_overfit.R` (synthetic sets, ties, thresholds, defaults; the
    16 DRC indicators against the saved set B selection, using `tests/testthat/fixtures/`).

* **`xgb()` truncates the lower interval bounds of a log-transformed outcome at 0.** Intervals
  were truncated to [0, 1] for `arcsin` and at 0 for `poisson` and gamma outcomes, but not for
  `transformation = "log"`, so on the natural interval scale `Lower` and `Lower_bench` of a
  positive outcome such as the Theil index could be negative. They are now `pmax(., 0)`, as for
  `poisson`. Only bounds that were negative change (to 0); everything else, including
  `interval_scale = "transformed"` (positive by construction), is identical. `xgb()` already stops on a
  non-positive outcome for `"log"` (`xgb_check1()`), so a log outcome is always positive and the truncation
  is always valid; `"log.shift"` is untouched. `sqrt` is not
  clipped (an open question, not changed here).
  - Test: `tests/testthat/test_xgb_log_truncation.R`.

* **`xgb()` no longer lets the benchmark-target perturbation blow up for a group pinned at the boundary (arcsine).**
  With `perturb_benchmark = TRUE` the rate-scale standard error of a benchmark target is turned
  into an increment on the transformed scale by the delta-method slope dz/dp, evaluated at the
  group's mean benchmarked rate. For the arcsine, dz/dp = 1 / (2 sqrt(p (1 - p))) is unbounded
  at 0 and 1, and the rate was floored at 1e-9, giving a slope of about 100 (a finite
  difference over 1e-4) for a group benchmarked to a target of exactly 0. With a zero standard
  error that was harmless; with the positive, boundary-corrected one the guard now asks for
  (Agresti-Coull), the perturbation became a fraction of a radian to a radian: about half of
  the replicates ended at 0, the rest spread over (0, 1), and the upper bounds of the group's
  domains reached 0.1 to 0.6 (found in the DRC audit, ROUND3 section 4). The slope is now
  evaluated at no less than the target standard error from either bound: the rate is clamped
  into [se, 1 - se] (arcsine only; at se of 0.5 or more the slope is taken at 0.5). **Results
  change only for a group whose mean benchmarked rate lies within one standard error of 0 or 1
  and whose standard error is positive**; for every other group, and when the standard error is
  0 or missing, the slope, and so the intervals, are identical to before. The log and square-root
  transformations are unchanged: their slopes (1/p, 1/(2 sqrt(p))) are unbounded at 0 too, but
  only a log-scale outcome with a mean near 0 and a positive standard error would be affected,
  and none of the indicators estimated so far is.
  - Test: `tests/testthat/test_xgb_perturb_slope.R` (the helper is unchanged away from the
    boundary; through `xgb()` the realised slope equals the old one for groups away from the
    boundary and the clamped one, about 7 rather than about 100, for a group pinned at zero,
    whose upper bounds stay of the order of the standard error).

* **`xgb()` builds the benchmarked interval with the basic-bootstrap sign by default (`interval_method = "basic"`).**
  A replicate's error is its benchmarked estimate minus its simulated truth, centred. The
  interval was built as the estimate *plus* the quantiles of that error, i.e. [estimate + q(a/2),
  estimate + q(1 - a/2)]. For an error defined as estimate minus truth the basic bootstrap
  interval for the truth *subtracts* the reversed quantiles, [estimate - q(1 - a/2), estimate -
  q(a/2)], the correction megb received earlier in this release. The two agree when the errors are
  symmetric about zero and differ when they are skewed, typically near a bound. `"basic"` is now the
  default; `interval_method = "plus"` reproduces the results of earlier versions exactly. It applies on either `interval_scale`.
  The unbenchmarked interval is the spread of the simulated truth around the estimate, not an
  error interval, and is not affected. Simulation (100 datasets at moderate rates, 60 at low
  rates): benchmarked coverage 0.707 -> 0.709 (moderate), 0.778 -> 0.790 on the natural scale and
  0.804 -> 0.812 on the transformed scale (low rates).
  - Test: `tests/testthat/test_xgb_interval_method.R`.

* **`xgb()` stops on a zero, negative or non-finite `benchmark_target_se` when
  `perturb_benchmark = TRUE`.** A direct estimate on the boundary of its range (a proportion of
  exactly 0 or 1, a sample whose clusters all agree) has a design-based standard error of exactly
  zero. The perturbation of the benchmark target then added nothing for that group, and the
  intervals of its domains collapsed towards the target, silently. The error asks for a
  boundary-corrected standard error (for a proportion at 0 or 1, Agresti-Coull, which needs the
  effective sample size that only the caller has); the package does not substitute one. A zero
  can be accepted deliberately with the new `allow_zero_benchmark_se = TRUE` (a warning names the
  groups); negative and non-finite values always stop. The internal Horvitz-Thompson standard
  error used when no `benchmark_target_se` is supplied is checked the same way (it is zero when
  every sampled unit of a benchmark group has the same outcome, or the group has fewer than two
  units). With `perturb_benchmark = FALSE` no standard error is used and nothing is checked.
  - Test: `tests/testthat/test_xgb_benchmark_se_guard.R`.

* **`xgb()` can build its intervals on the transformation scale (`interval_scale = "transformed"`).**
  The default (`"natural"`) is unchanged: quantiles of the bootstrap error on the scale of the
  estimates, added to the estimate, then truncated to [0, 1] for `arcsin` (and at 0 for `poisson`).
  Near a bound that piles intervals onto 0 or 1 (with the benchmark-target perturbation, half the
  wards of a poverty map had an upper bound of exactly 1), and for `log` outcomes such as an
  inequality index it gives negative lower bounds. With `"transformed"` the benchmarked replicates
  and the simulated truths are transformed, the quantiles of their centred differences are added to
  the transformed benchmarked estimate and the result is back-transformed; the unbenchmarked
  interval is the same construction on the point-estimate replicates. Bounds are valid by
  construction and nothing is clipped. Point estimates and variances do not change. Available for
  `transformation` `"no"`, `"arcsin"`, `"log"` and `"sqrt"`.
  - New opt-in `keep_replicates = TRUE` stores the per-replicate domain estimates (point,
    benchmarked before and after the benchmark-target perturbation), the simulated truths and the
    perturbations in `$replicates`, so intervals can be rebuilt on any scale without re-running the
    bootstrap. Off by default (B x domains entries per matrix).
  - Test: `tests/testthat/test_xgb_interval_scale.R`.

* **megb's bootstrap is centred on the cross-fitted fit.** With
  `predict_sampled = "crossfit"`, each sampled domain is predicted by the fold model
  that held it out, and so is the bootstrap truth. The bootstrap outcome, however, was
  centred on the full booster's in-sample fit, and its residual pool was taken around
  that fit plus the random effect. Both now use the fold models' out-of-fold
  predictions instead. The random effect is the point fit's, as before. This applies
  to every `bootstrap_refit` option, including the default `"leaves_only"`. With
  `predict_sampled = "full"` nothing changes.
  - The full-refit bootstrap's replicate fits now use the point fit's EM settings
    (25 iterations, tolerance 1e-4), defined in one place
    (`.megb_em_max_iterations`). They were capped at 10.
  - Point estimates are unchanged. Intervals change under `predict_sampled = "crossfit"`
    (the default since 1.0.1), and under `bootstrap_refit = "full"`.
  - Test: `tests/testthat/test-megb-bootstrap-centre.R`.

* **megb's interval had its bootstrap error quantiles on the wrong side.** The interval
  is built from each replicate's prediction error, estimate minus bootstrap truth. It
  was the point estimate *plus* the error quantiles. The basic bootstrap interval
  *subtracts* them: [estimate - q(1 - alpha/2), estimate - q(alpha/2)]. The two agree only
  when the errors are symmetric. With skewed errors the old interval put its long tail
  on the wrong side of the estimate. This applies to every `mse_type = "var"` interval,
  benchmarked or not. The variance (of the centred errors) is unchanged.

* **The `leaves_only` bootstrap can run the point fit's EM in each replicate (opt-in).**
  By default each replicate still refreshes the leaves once on the replicate outcome,
  random effect included, and fits the random effect once. This single pass is a
  documented approximation to the point fit's EM.
  - `options(povmap.leaves_only.em = TRUE)` runs the EM with the trees fixed. It
    starts from a random-effect-only fit on the replicate outcome, then alternates:
    the fold models' leaves are refreshed on the outcome net of the estimated random
    effect, and the random effect is re-estimated by ML from those models'
    out-of-fold residuals. It stops at the point fit's iteration cap and tolerance.
  - `options(povmap.leaves_only.em_start = "point")` starts the EM from the point
    fit's random effect instead.
  - Why the single pass is the default: on the Colombian M19 gate (six outcomes,
    B = 100, same seeds), the EM took a median of 4 iterations (range 3 to 6) and
    cost 3.3 times the single pass (2.8 times with the warm start). It changed
    coverage by -0.4 points in sample and -0.2 out of sample, and mean width by
    -0.1 and -0.7 percent.
  - `boot_diag$boot_em_iterations` reports each replicate's iteration count (0 for
    the single pass).
  - Test: `tests/testthat/test-megb-bootstrap-interval-em.R`.

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