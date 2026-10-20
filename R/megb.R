#' Mixed Effects Gradient Boosting for domain-level averages
#'
#' The function \code{megb} combines gradient boosting with linear mixed models
#' (via an EM algorithm) to estimate domain-level averages for small area
#' estimation. It wraps the \code{megb} function from the \pkg{MEGB} package,
#' returning a \code{povmap} object compatible with \code{write.excel},
#' \code{estimators}, and \code{map_plot}.
#'
#' Unlike \code{\link{xgb}}, MEGB operates directly at the unit level and
#' models area-level random effects explicitly through a linear mixed model.
#' No sub-domain identifier is required: the EM algorithm handles the
#' sample-to-population matching internally.
#'
#' @param fixed a two-sided linear formula describing the model. The dependent
#'   variable is on the left of \code{~} and the predictors on the right,
#'   separated by \code{+}. The domain identifier named in \code{domains} may
#'   appear in the formula (it is automatically excluded from the GB predictor
#'   matrix and used only as the LMM grouping variable).
#' @param smp_data a data frame containing all variables in \code{fixed} plus
#'   the domain identifier.
#' @param smp_weights character string naming the survey-weight variable in
#'   \code{smp_data}. Rescaled within each domain to average 1, the weights
#'   weight the booster's training rows, the linear mixed model fits, and the
#'   bootstrap's leaf refresh and residual draws. Defaults to \code{NULL}.
#' @param pop_data a data frame containing the domain identifier and all
#'   predictor variables from \code{fixed}.
#' @param pop_weights character string naming a population-count or expansion-
#'   weight variable in \code{pop_data}. When supplied, domain means are
#'   computed as weighted averages of unit-level predictions. Defaults to
#'   \code{NULL} (simple mean).
#' @param domains character string naming the domain variable in both
#'   \code{smp_data} and \code{pop_data}.
#' @param transformation transformation applied to the dependent variable before
#'   fitting. One of \code{"no"} (default), \code{"log"}, \code{"log.shift"},
#'   \code{"arcsin"}, \code{"logistic"}, or \code{"poisson"}. Unit-level
#'   population predictions are back-transformed before domain aggregation.
#' @param mse logical. If \code{TRUE} (default), runs the parametric bootstrap
#'   to quantify uncertainty. The form of the resulting summary is controlled
#'   by \code{mse_type}.
#'   The EM loop runs to convergence (relative log-likelihood change <
#'   1e-04) or for at most 25 iterations, whichever comes first. If
#'   \eqn{\sigma_u} collapses toward zero, that typically reflects
#'   identifiability between area-level covariates (which are constant
#'   within an area, and so collinear with the area indicator) and the
#'   area random effect, rather than a failure of EM. Diagnose by
#'   comparing specifications with and without the area-level covariates
#'   and by leave-areas-out cross-validation, rather than by truncating
#'   the EM loop.
#' @param bench_target one of \code{"fixed"} (default), \code{"random"} or \code{"truth"}.
#'   Controls how the benchmark target is treated *inside the bootstrap*.
#'   \describe{
#'     \item{\code{"fixed"}}{rakes every replicate to the observed (province/level)
#'       target, so the benchmarked MSE/CI is the spread of the estimate around the
#'       realised truth at the observed target. This is the default and the choice
#'       for benchmarking to an observed survey target.}
#'     \item{\code{"random"}}{recomputes the target from each replicate's own resampled
#'       survey, weighted by \code{benchmark_weights}. The benchmarked estimate then
#'       tracks its own bootstrap mean and the level variance collapses, giving
#'       too-narrow intervals. Only matters when \code{benchmark} is character
#'       (internal, survey-derived benchmarking); for numeric or data.frame
#'       benchmarks the target is a constant by construction.}
#'     \item{\code{"truth"}}{rakes replicate b to \eqn{T^*_{g,b} = \tau_{g,b} + \epsilon_{g,b}},
#'       where \eqn{\tau_{g,b}} is the simulated truth of benchmark group g in that
#'       replicate, aggregated from the domain truths with the population totals the
#'       benchmark aggregates the estimate with, and \eqn{\epsilon_{g,b} \sim N(0, SE_g^2)}
#'       is the error of the observed target, with \eqn{SE_g} from \code{benchmark_target_se}.
#'       The provincial error of the benchmarked estimate is then \eqn{\epsilon_{g,b}}
#'       alone: the sampling error of the target is counted once, and the spread
#'       of the simulated provincial truth is not added to it. The \eqn{B \times G}{B x G}
#'       draws come once from \code{seed} on a stream of their own (the same stream and
#'       draws as \code{perturb_benchmark}), the target is clamped like the perturbed
#'       replicates (lower bound 1e-9 for arcsin, log and poisson, upper bound 1 - 1e-9
#'       for arcsin), and the point estimate is unchanged, raked to the observed
#'       \code{benchmark}. With \eqn{SE_g = 0} the provincial error is zero. Needs
#'       \code{benchmark}, \code{benchmark_level}, \code{pop_weights}, \code{mse = TRUE}
#'       and \code{benchmark_target_se}, and excludes \code{perturb_benchmark}.
#'       Works with every \code{bootstrap_refit}. The result then holds
#'       \code{bench_target = "truth"} and \code{benchmark_target_draws}.}
#'   }
#' @param em_maxit the maximum number of EM iterations of the point fit (and of each
#'   replicate fit of \code{bootstrap_refit = "full"}). Default 25. The EM stops earlier
#'   when the relative log-likelihood change falls below \code{em_tol}. A resampled
#'   cross-validation inside each iteration makes the log-likelihood jitter, so with
#'   boosting the rule is often not met and the fit runs to \code{em_maxit}; check the
#'   \eqn{\sigma_b^2} trajectory (\code{megb_model$iterations_used}) before relying on it.
#' @param em_tol the EM stopping tolerance on the relative change of the log-likelihood.
#'   Default \code{1e-04}.
#' @param perturb_benchmark logical, default \code{FALSE}. If \code{TRUE}, the benchmarked
#'   bootstrap replicates are perturbed by the sampling error of the benchmark target,
#'   exactly as in \code{\link{xgb}}: after the bootstrap, a draw from
#'   \eqn{N(0, SE_g^2)} is added to the benchmarked replicate of every domain in benchmark
#'   group \eqn{g}, on the model's transformed scale (the increment is the draw times the
#'   local slope of the transformation, \code{.perturb_slope}) and back-transformed; for
#'   \code{transformation = "no"} it is added directly. The truth side is not perturbed. The
#'   \eqn{B \times G}{B x G} draws are made once from \code{seed} (a stream of their own, so the
#'   bootstrap is unchanged). Needs \code{benchmark}, \code{benchmark_level},
#'   \code{mse = TRUE} and \code{benchmark_target_se}. With \code{FALSE} the result is
#'   identical to povmap 1.0.3. The result then also holds \code{Var_bench_unperturbed} and
#'   \code{benchmark_perturbation}. Exclusive with \code{bench_target = "truth"}, which draws
#'   the target error itself.
#' @param benchmark_target_se a named numeric vector of the standard errors of the
#'   benchmark targets, keyed by \code{benchmark_level} group; the SD of the perturbation
#'   (\code{perturb_benchmark = TRUE}) or of the target error (\code{bench_target = "truth"}).
#'   Required by either; a missing group, or a non-finite,
#'   negative or (unless \code{allow_zero_benchmark_se}) zero value is an error. Unlike
#'   \code{xgb}, there is no internal Horvitz-Thompson fallback.
#' @param allow_zero_benchmark_se logical, default \code{FALSE}. Allows a zero standard
#'   error for some benchmark group (a warning is given; the perturbation or the target
#'   error adds nothing there).
#' @param arcsin_backtransform one of \code{"exact"} (default) or \code{"plain"}; used only with
#'   \code{transformation = "arcsin"}. \code{"exact"} back-transforms each census cell with the exact
#'   mean of \eqn{\sin^2} under normal error, \eqn{1/2 - 1/2 \cos(2\mu) e^{-2v}}, where \eqn{v} is
#'   \eqn{\sigma_e^2} plus the conditional variance of the area random effect. \code{"plain"} uses
#'   \eqn{\sin^2(\mu)} with no variance term, everywhere the back-transformation is used: the point
#'   estimate, the bootstrap area estimates and truths, the benchmark target of the replicates, and the
#'   random-effect variance term (which is then not computed). REStat's preregistered census-truth test
#'   found that the exact form overshoots near zero and biases low-prevalence areas upward; the plain
#'   form raised sampled-area coverage by 0.06 to 0.16 and cut that bias by about three quarters at no
#'   RMSE cost. The default is unchanged for now. The result then holds \code{arcsin_backtransform =
#'   "plain"}.
#' @param boot_synthetic logical, default \code{FALSE}. \code{TRUE} (needs \code{mse = TRUE} and \code{bootstrap_refit = "leaves_only"}) adds
#'   \code{boot_diag$synthetic}: the point fit's cross-fitted synthetic area prediction (booster part only, no area effect, the back-transformation of an
#'   unsampled area, unbenchmarked, on the reporting scale), the same prediction of every bootstrap replicate, the replicates' unbenchmarked truth and the
#'   bootstrap MSE of the synthetic prediction per area. Estimates, MSE and intervals are unchanged (bit for bit). Input of \code{megb_level_calibration()}.
#' @param boot_area_var how the bootstrap world's area-level variance is built; one of
#'   \code{"legacy"} (default), \code{"centred"} or \code{"centred_mm"}. Point estimates are the
#'   same under every option; only the bootstrap (hence the MSE and the intervals) changes.
#'   \code{"legacy"} reproduces the previous behaviour bit for bit: each sampled area redraws its own
#'   residuals (rescaled globally to the fit's error SD), and the bootstrap area effects are drawn from
#'   the rescaled EBLUPs. \code{"centred"}: the residual pool of each area is centred within the area
#'   (mean over the area's sampled units, the weighted mean when \code{weightedBS = TRUE}) and scaled by
#'   \eqn{\sqrt{n_i / (n_i - 1)}}; an area with a single sampled unit, and the unsampled areas, draw
#'   from the pooled centred, scaled residuals of the areas with at least two units (single-unit areas
#'   excluded). \code{"centred_mm"}: \code{"centred"} plus bootstrap area effects drawn, with
#'   replacement, from the area deviations \eqn{d_i} (the weighted mean of the response minus its
#'   cross-fitted prediction in each sampled area, centred with the area weights) rescaled to the
#'   standard deviation \eqn{\hat\sigma_b(MM)}, where \eqn{\hat\sigma_b^2(MM) = \max(0, var_w(d_i) -
#'   mean_w(\hat\sigma_e^2 / n_i))} with area weights the sums of the unit weights and, when weights are used, the effective area size (sum w)^2 / sum w^2 in the noise term; zero gives no area
#'   effect. See the preregistration MEGB_BOOT_FIX_PREREG_v4.md.
#' @param bootstrap_refit one of \code{"leaves_only"} (default),
#'   \code{"full"}, or \code{"lmm_only"}.
#'   \code{"leaves_only"} (xgboost only) freezes the tree structure of the
#'   original GB across all bootstrap iterations and only refreshes leaf-node
#'   values via xgboost's \code{refresh} updater. This captures within-model
#'   leaf-prediction sampling uncertainty plus the LMM random-effect
#'   uncertainty, without exposing the SE to GB tree-selection variance —
#'   which empirically overstates uncertainty for sparse or boundary-near
#'   indicators (e.g. ownership rates, livestock counts, access indicators
#'   with many near-zero territories) where bootstrap-perturbed residuals
#'   drive fresh GB fits to substantially different trees. For dense
#'   indicators (continuous welfare, well-spread proportions) \code{"leaves_only"}
#'   gives results within ~20\% of \code{"full"}; for sparse indicators the
#'   two can differ by an order of magnitude. Recommended default and aligned
#'   with how most SAE methods report SEs (conditional on chosen model
#'   structure). Each replicate refreshes the leaves once on the replicate
#'   outcome and fits the random effect once, an approximation to the point
#'   fit's EM; \code{options(povmap.leaves_only.em = TRUE)} runs that EM with
#'   the trees fixed (about three times the cost; see NEWS).
#'   \code{"full"}: refit both the gradient booster and the linear mixed model
#'   in every bootstrap iteration. This additionally captures GB
#'   tree-selection variance and is appropriate when that is the inferential
#'   target (e.g. inference about a model \emph{class} rather than a fitted
#'   model). For sparse indicators this can substantially overstate true
#'   predictive uncertainty because the residual bootstrap drives fresh GBs
#'   into pathological refits.
#'   \code{"lmm_only"}: refit only the linear mixed model on bootstrap residuals
#'   each iteration, treating the gradient-booster fit as fixed. This is the
#'   classical EBLUP parametric bootstrap (Prasad–Rao / Hall–Maiti) and yields
#'   the textbook leading-order MSE \eqn{g_1 = \gamma_d \sigma_e^2 / n_d}.
#'   It runs faster than the other modes but assumes the FE estimator is
#'   parametric and contributes negligible variance — an assumption that
#'   holds for low-dimensional linear FE models but \emph{not} for boosted
#'   trees, where it typically understates SE by an order of magnitude. Use
#'   only if you explicitly want the classical-EBLUP variance decomposition.
#' @param mse_type one of \code{"var"} (default) or \code{"mse"}. \code{"var"}
#'   reports the variance of the bootstrap prediction errors (each replicate's
#'   domain mean minus its bootstrap truth, centred per domain) and forms the CI
#'   as the basic bootstrap interval: the point estimate minus the errors' upper
#'   and lower quantiles (the quantiles were added before 5 Oct 2026). Benchmarked CIs use the centred residual
#'   between benchmarked-bootstrap and bootstrap-truth, again as in \code{xgb}.
#'   \code{"mse"} reports the prediction MSE \eqn{E[(\hat\theta-\theta)^2]} —
#'   computed in transformed space and back-mapped via the delta method for the
#'   unbenchmarked estimator, and directly in original space for the
#'   benchmarked estimator — and uses a normal-approximation CI. \code{"mse"}
#'   is generally wider than \code{"var"} because it includes the irreducible
#'   variability of the random effect.
#' @param B number of parametric bootstrap iterations for MSE estimation.
#'   Defaults to \code{100}.
#' @param bootstrap_cores number of cores for parallel bootstrap. \code{0} or
#'   \code{1} runs sequentially. Defaults to \code{0}.
#' @param conf_level confidence level for the confidence interval. Defaults to
#'   \code{0.95}.
#' @param gradient_params named list of gradient boosting hyperparameters passed
#'   to the chosen \code{gbm_engine}. When \code{NULL}, engine defaults are
#'   used. For XGBoost the relevant names are \code{eta}, \code{nrounds},
#'   \code{max_depth}, \code{subsample}, etc.
#' @param gbm_engine gradient boosting engine: \code{"xgboost"} (default),
#'   \code{"lightgbm"}, or \code{"catboost"}.
#' @param na.rm if \code{TRUE}, rows with \code{NA} values are removed from
#'   both \code{smp_data} and \code{pop_data} before fitting. Defaults to
#'   \code{FALSE}.
#' @param seed integer seed for reproducibility. Defaults to \code{123}.
#' @param cv_folds how the internal cross-validation (\code{xgboost::xgb.cv},
#'   xgboost engine only) splits the sample in every EM iteration. Its
#'   out-of-fold predictions give the residuals from which the linear mixed
#'   model estimates the error variance \eqn{\sigma_e^2}, the random-effect
#'   variance \eqn{\sigma_u^2} and the random effects, and the fold models set
#'   the number of boosting rounds when \code{early_stopping = TRUE}.
#'   \code{"domain"} (default) holds out \code{cv_nfold} folds of whole domains, so no domain
#'   has rows on both sides of a fold. \code{"rows"} holds out \code{cv_nfold} folds of random
#'   rows, the behaviour before povmap 6bd056f.
#'
#'   Domain folds matter most when covariates are constant within domains (for
#'   example, municipal or area-level covariates attached to every household in
#'   the municipality). With row folds, each held-out row's domain is
#'   represented in the training folds by its other rows with identical
#'   covariate values, so the booster can learn the domain's mean outcome and the
#'   out-of-fold residuals carry little between-domain variation. The random-
#'   effect variance is then estimated at or near zero, the random effects
#'   vanish, and the bootstrap intervals, whose out-of-sample width is driven by
#'   \eqn{\sigma_u}, become far too narrow. Holding out whole domains makes
#'   the out-of-fold residuals honest about between-domain variation that the
#'   covariates do not explain, which is the variation the random effect is
#'   meant to capture. In the Colombia migration evaluation, row folds gave
#'   \eqn{\hat\sigma_u = 0} in about half the fits on a covariate set in which
#'   half the covariates are municipality means, with out-of-sample coverage of
#'   95\% intervals of 0.44 to 0.53; domain folds gave about 0.88 at the same
#'   point accuracy.
#'
#'   Domain labels come from \code{domains}. If \code{"domain"} is requested
#'   but no labels reach the booster, it falls back to row folds with a
#'   message. When not supplied, the default is taken from
#'   \code{getOption("megb.cv_folds", "domain")}.
#' @param early_stopping logical. If \code{TRUE} (default), the internal
#'   cross-validation stops when the held-out error has not improved for 10
#'   rounds, and the full-data booster is trained to that round. If
#'   \code{FALSE}, the fold models and the booster train to the full
#'   \code{nrounds} in \code{gradient_params}, as \code{\link{xgb}} does.
#'   When not supplied, the default is taken from
#'   \code{getOption("megb.early_stopping", TRUE)}.
#' @param predict_sampled how the booster part of a sampled domain's prediction
#'   is formed. \code{"crossfit"} uses the fold model of the final EM iteration
#'   that held the domain out; \code{"full"} uses the booster trained on the
#'   whole sample. The random effects are estimated from the fold models'
#'   out-of-fold residuals, so under \code{"crossfit"} the booster part and the
#'   random effect are measured against the same baseline. Under \code{"full"}
#'   the booster has seen the domain's own rows, which, with covariates constant
#'   within domains (area-level covariates), can count part of the domain's
#'   deviation twice. The default, \code{NULL}, is \code{"crossfit"} wherever it
#'   applies (xgboost, \code{cv_folds = "domain"}, and with \code{mse = TRUE} a
#'   \code{"leaves_only"} or \code{"full"} bootstrap) and \code{"full"}, with a
#'   message, otherwise; asking for \code{"crossfit"} where it does not apply is
#'   an error. Until povmap 1.0.1 the default was \code{"full"}. In a Colombian
#'   evaluation with 19 and 24 municipal covariates, cross-fitting with ten folds
#'   and \code{predict_unsampled = "foldmean"} lowered RMSE by 4.3 and 3.5
#'   percent against \code{"full"} with five folds (8.4 and 6.4 percent in
#'   sample); with sub-area covariates it was neutral. The fit itself, including
#'   the random effects and variance components, does not depend on this
#'   argument. The result's \code{$crossfit$ind} gives the estimates of the same
#'   fit under the other settings: \code{Mean_full} (full booster everywhere),
#'   \code{Mean_crossfit} (cross-fitted sampled, full unsampled) and
#'   \code{Mean_foldmean} (cross-fitted sampled, fold-model average unsampled).
#'   With \code{mse = TRUE} the bootstrap predicts its replicates the same way:
#'   the \code{"leaves_only"} bootstrap from each refreshed fold model, the
#'   \code{"full"} bootstrap from each replicate fit's own fold models.
#' @param predict_unsampled how the booster part of an unsampled domain's
#'   prediction is formed. \code{"foldmean"} averages the predictions of the fold
#'   models of the final EM iteration; \code{"full"} uses the booster trained on
#'   the whole sample. Unsampled domains have no random effect either way. The
#'   default, \code{NULL}, is \code{"foldmean"} wherever fold-model prediction
#'   applies (as for \code{predict_sampled}) and \code{"full"}, with a message,
#'   otherwise. Until povmap 1.0.1 unsampled domains always used the full
#'   booster. In the Colombian evaluation the fold-model average lowered
#'   out-of-sample RMSE by 0.3 to 0.8 percent. With \code{mse = TRUE} the
#'   \code{"leaves_only"} bootstrap averages the refreshed fold models, and the
#'   \code{"full"} bootstrap each replicate fit's fold models.
#' @param cv_nfold the number of folds of the internal cross-validation
#'   (xgboost only), of whole domains or of rows according to \code{cv_folds}.
#'   Default 10 (5 until povmap 1.0.1). More folds train each fold model on a
#'   larger share of the sample: with fold-model prediction, ten folds beat five
#'   by 0.3 to 0.5 percent in sample in the Colombian evaluation. Each EM
#'   iteration costs proportionally more: a point fit took 1.6 to 1.9 times as
#'   long with ten folds as with five, and a fit with the default bootstrap
#'   (B = 100) about 2.1 times as long as the five-fold, full-prediction
#'   configuration. Under fold-model prediction there is one fold model per fold.
#'   When not supplied, the default is taken from
#'   \code{getOption("megb.cv_nfold", 10L)}. For the settings before povmap
#'   1.0.1 use \code{cv_nfold = 5, predict_sampled = "full",
#'   predict_unsampled = "full"}.
#' @param ... additional arguments forwarded to \code{MEGB::megb}.
#'
#' @return An object of class \code{c("megb", "xgb", "povmap")} with elements:
#'   \describe{
#'     \item{\code{ind}}{data frame of domain-level point estimates
#'       (\code{Domain}, \code{Mean}).}
#'     \item{\code{var}}{data frame of uncertainty estimates (\code{Domain},
#'       \code{Mean}): with \code{mse_type = "var"} the variance of the centred
#'       bootstrap prediction errors, with \code{"mse"} the delta-method-corrected
#'       MSE; \code{NULL} if \code{mse = FALSE}.}
#'     \item{\code{CI}}{data frame of confidence intervals (\code{Domain},
#'       \code{Lower}, \code{Upper}): with \code{mse_type = "var"} the basic
#'       bootstrap interval from the prediction-error quantiles, with \code{"mse"}
#'       a normal approximation; \code{NULL} if \code{mse = FALSE}.}
#'     \item{\code{yhat}}{data frame of back-transformed sample-unit predictions
#'       (\code{obs_id}, \code{hat}).}
#'     \item{\code{model}}{the fitted GB model object (xgboost booster when
#'       \code{gbm_engine = "xgboost"}).}
#'     \item{\code{megb_model}}{the full MEGB model list including both the
#'       booster and the LME4 random-effects model.}
#'     \item{\code{gbm_engine}}{the engine used.}
#'     \item{\code{smp_data}}{the original sample data.}
#'     \item{\code{out_call}}{the matched call.}
#'     \item{\code{transformation}}{the transformation string.}
#'     \item{\code{framework}}{internal framework list with domain metadata.}
#'   }
#'
#' @seealso \code{\link{xgb}}, \code{\link{write.excel}},
#'   \code{\link{estimators}}
#'
#' @references
#' Merfeld, J. D., Dang, H., & Newhouse, D. (2025). Improving Estimates of
#' Mean Welfare and Uncertainty in Developing Countries. The World Bank.
#'
#' @importFrom stats qnorm tapply
#' @export
#'
#' @examples
#' \donttest{
#' data("eusilcA_pop")
#' data("eusilcA_smp")
#'
#' megb_model <- megb(
#'   fixed   = eqIncome ~ eqsize + cash + self_empl + unempl_ben +
#'             age_ben + surv_ben + sick_ben + dis_ben + rent +
#'             fam_allow + house_allow + cap_inv + tax_adj,
#'   smp_data = eusilcA_smp,
#'   pop_data = eusilcA_pop,
#'   domains  = "district",
#'   gradient_params = list(eta = 0.1, nrounds = 100, max_depth = 3)
#' )
#'
#' estimators(megb_model, indicator = "Mean", MSE = TRUE, CV = TRUE)
#' }

megb <- function(fixed,
                 smp_data,
                 smp_weights      = NULL,
                 pop_data,
                 pop_weights      = NULL,
                 domains,
                 transformation   = "no",
                 mse              = TRUE,
                 mse_type         = c("var", "mse"),
                 bootstrap_refit  = c("leaves_only", "full", "lmm_only"),
                 bench_target     = c("fixed", "random", "truth"),
                 B                = 100,
                 bootstrap_cores  = 0,
                 conf_level       = 0.95,
                 gradient_params  = NULL,
                 gbm_engine       = "xgboost",
                 na.rm            = FALSE,
                 seed             = 123,
                 benchmark        = NULL,
                 benchmark_type   = "ratio",
                 benchmark_level  = NULL,
                 benchmark_weights = NULL,
                 weightedBS       = TRUE,
                 cpus             = NULL,
                 cv_folds         = getOption("megb.cv_folds", "domain"),
                 early_stopping   = getOption("megb.early_stopping", TRUE),
                 predict_sampled  = NULL,
                 predict_unsampled = NULL,
                 cv_nfold         = getOption("megb.cv_nfold", 10L),
                 em_maxit         = 25L,
                 em_tol           = 1e-04,
                 perturb_benchmark = FALSE,
                 benchmark_target_se = NULL,
                 allow_zero_benchmark_se = FALSE,
                 arcsin_backtransform = c("exact", "plain"),
                 boot_area_var = c("legacy", "centred", "centred_mm"),
                 boot_synthetic = FALSE,
                 ...) {
  # weightedBS: when TRUE and smp_weights is supplied, the bootstrap residual
  # sampling uses probabilities proportional to weights (mirroring xgb's
  # weightedBS semantics). When FALSE, residuals are sampled with equal
  # probability (the historical default before weights were plumbed through).

  out_call        <- match.call()
  mse_type        <- match.arg(mse_type)
  bootstrap_refit <- match.arg(bootstrap_refit)
  bench_target    <- match.arg(bench_target)
  arcsin_backtransform <- match.arg(arcsin_backtransform)
  boot_area_var   <- match.arg(boot_area_var)
  boot_synthetic  <- isTRUE(boot_synthetic)
  if (boot_synthetic && !mse) stop("boot_synthetic = TRUE needs mse = TRUE.")
  cv_folds        <- match.arg(cv_folds, c("domain", "rows"))
  if (!is.logical(early_stopping) || length(early_stopping) != 1L || is.na(early_stopping))
    stop("early_stopping must be TRUE or FALSE.")
  if (!is.numeric(cv_nfold) || length(cv_nfold) != 1L || is.na(cv_nfold) || cv_nfold < 2 || cv_nfold != round(cv_nfold))
    stop("cv_nfold must be a whole number of at least 2.")
  cv_nfold <- as.integer(cv_nfold)
  if (!is.numeric(em_maxit) || length(em_maxit) != 1L || is.na(em_maxit) || em_maxit < 1 || em_maxit != round(em_maxit))
    stop("em_maxit must be a whole number of at least 1.")
  em_maxit <- as.integer(em_maxit)
  if (!is.numeric(em_tol) || length(em_tol) != 1L || is.na(em_tol) || em_tol <= 0)
    stop("em_tol must be a positive number.")
  if (!is.logical(perturb_benchmark) || length(perturb_benchmark) != 1L || is.na(perturb_benchmark))
    stop("perturb_benchmark must be TRUE or FALSE.")
  if (perturb_benchmark) {
    if (is.null(benchmark) || is.null(benchmark_level) || !mse)
      stop("perturb_benchmark = TRUE needs a benchmark, a benchmark_level and mse = TRUE.")
    if (is.null(benchmark_target_se))
      stop("perturb_benchmark = TRUE needs benchmark_target_se, a named vector of the standard errors of the benchmark targets, keyed by benchmark_level group.")
  }
  if (bench_target == "truth") {
    if (perturb_benchmark)
      stop("bench_target = \"truth\" already draws the target error (it replaces perturb_benchmark); set perturb_benchmark = FALSE.")
    if (is.null(benchmark) || is.null(benchmark_level) || !mse)
      stop("bench_target = \"truth\" needs a benchmark, a benchmark_level and mse = TRUE.")
    if (is.null(benchmark_target_se))
      stop("bench_target = \"truth\" needs benchmark_target_se, a named vector of the standard errors of the benchmark targets, keyed by benchmark_level group.")
    if (is.null(pop_weights))
      stop("bench_target = \"truth\" needs pop_weights: the simulated truth of a benchmark group is aggregated with the weights the benchmark uses.")
  }
  # Fold-model prediction (predict_sampled = "crossfit", predict_unsampled =
  # "foldmean"; the defaults from povmap 1.0.1) needs the xgboost engine, folds of
  # whole domains and, with mse = TRUE, a bootstrap that refits it
  # ("leaves_only" or "full"). NULL means: the fold-model setting where it applies,
  # otherwise "full" with a message. An explicit request it cannot honour is an error.
  .fold_why <- c(if (gbm_engine != "xgboost") "gbm_engine is not \"xgboost\"",
                 if (cv_folds != "domain") "cv_folds is not \"domain\"",
                 if (mse && bootstrap_refit == "lmm_only") "bootstrap_refit = \"lmm_only\" has no fold-model bootstrap")
  .resolve <- function(x, arg, alt) {
    if (is.null(x)) {
      if (!length(.fold_why)) return(alt)
      message("megb: ", arg, " = \"full\" because ", paste(.fold_why, collapse = " and "), ".")
      return("full")
    }
    x <- match.arg(x, c("full", alt))
    if (x == alt && length(.fold_why))
      stop(arg, " = \"", alt, "\" is not available: ", paste(.fold_why, collapse = " and "), ".")
    x
  }
  predict_sampled   <- .resolve(predict_sampled, "predict_sampled", "crossfit")
  predict_unsampled <- .resolve(predict_unsampled, "predict_unsampled", "foldmean")
  .fold_pred <- predict_sampled == "crossfit" || predict_unsampled == "foldmean"

  # Local %||% so the diagnostic message below isn't fragile to NULLs.
  `%||%` <- function(a, b) if (is.null(a)) b else a

  if (!is.null(cpus)) bootstrap_cores <- cpus

  if (is.null(benchmark_weights) && !is.null(smp_weights))
    benchmark_weights <- smp_weights

  # ── 1. Framework ─────────────────────────────────────────────────────────────
  fwk <- framework_megb(
    fixed             = fixed,
    smp_data          = smp_data,
    pop_data          = pop_data,
    smp_weights       = smp_weights,
    pop_weights       = pop_weights,
    domains           = domains,
    transformation    = transformation,
    conf_level        = conf_level,
    na.rm             = na.rm,
    benchmark_weights = benchmark_weights
  )

  # ── 1b. Local benchmark helper ────────────────────────────────────────────────
  # Route through benchmark_xgb_level for level-based benchmarking so megb
  # supports the same benchmark_type set as xgb (notably "logit_raking",
  # which benchmark_ebp_level does not implement and would silently return
  # all-NA Mean_bench).
  add_benchmark_megb <- function(x, benchmark_level, fwk, fixed,
                                 benchmark, benchmark_type, domain_labels = NULL) {
    # Carry the domain labels so benchmark_xgb_level matches the estimates to domains by name (it refuses unlabelled
    # estimates). Every caller passes them: the point estimate and each bootstrap replicate (benchmark_fn).
    point_estim       <- if (!is.null(domain_labels))
                           list(ind = data.frame(Mean = x, Domain = domain_labels,
                                                 stringsAsFactors = FALSE))
                         else list(ind = data.frame(Mean = x))
    if (is.null(benchmark_level)) {
      point_estim$ind <- benchmark_ebp_national(
        point_estim    = point_estim,
        framework      = fwk,
        fixed          = fixed,
        benchmark      = benchmark,
        benchmark_type = benchmark_type)
    } else {
      point_estim$ind <- benchmark_xgb_level(
        point_estim     = point_estim,
        framework       = fwk,
        fixed           = fixed,
        benchmark       = benchmark,
        benchmark_type  = benchmark_type,
        benchmark_level = benchmark_level)
    }
    point_estim$ind$Mean_bench
  }
  # Lean env: only package functions needed, no local-scope capture
  environment(add_benchmark_megb) <- getNamespace("povmap")

  # ── 2. Transformation functions ───────────────────────────────────────────────
  if (transformation == "arcsin") {
    transform_outcome      <- arcsin_transform
    back_transform_outcome <- arcsin_transform_back
  } else if (transformation == "log") {
    transform_outcome      <- log_transform
    back_transform_outcome <- log_transform_back
  } else if (transformation == "log.shift") {
    ls_lambda              <- if (min(fwk$Y_smp) <= 0) abs(min(fwk$Y_smp)) + 1 else 0
    transform_outcome      <- function(y) list(y = log(y + ls_lambda), shift = NULL)
    back_transform_outcome <- function(y) exp(y) - ls_lambda
  } else if (transformation == "logistic") {
    transform_outcome      <- logit_transform_epsilon
    back_transform_outcome <- logit_transform_back
  } else if (transformation == "poisson") {
    if (min(fwk$Y_smp) < 0)
      stop("Outcome must be non-negative when using the poisson transformation")
    transform_outcome      <- poisson_transform
    back_transform_outcome <- poisson_transform_back
  } else {
    transform_outcome      <- no_transform
    back_transform_outcome <- no_transform_back
  }

  # ── 3. Transform outcome ──────────────────────────────────────────────────────
  Y_smp_t <- transform_outcome(fwk$Y_smp)$y

  # ── 4. Prepare data for MEGB::megb ───────────────────────────────────────────
  # Domain must be a factor; align factor levels across sample and population
  all_levels <- union(
    as.character(unique(fwk$smp_data[[domains]])),
    as.character(unique(fwk$pop_data[[domains]]))
  )

  smp_megb           <- fwk$smp_data[, c(domains, fwk$covariates), drop = FALSE]
  smp_megb[[domains]] <- factor(smp_megb[[domains]], levels = all_levels)

  pop_megb           <- fwk$pop_data[, c(domains, fwk$covariates), drop = FALSE]
  pop_megb[[domains]] <- factor(pop_megb[[domains]], levels = all_levels)

  X_smp <- smp_megb[, fwk$covariates, drop = FALSE]

  # ── 5. Fit MEGB model (MSE computed separately below after corrected_bt is known)
  megb_fit <- megb_em(
    Y               = Y_smp_t,
    X               = X_smp,
    dom_name        = domains,
    smp_data        = smp_megb,
    pop_data        = pop_megb,
    gradient_params = gradient_params,
    na.rm           = FALSE,
    seed            = seed,
    mse             = FALSE,
    gbm_engine      = gbm_engine,
    smp_weights_vec = fwk$smp_weights_vec,
    cv_folds        = cv_folds,
    early_stopping  = early_stopping,
    predict_sampled = predict_sampled,
    predict_unsampled = predict_unsampled,
    cv_nfold        = cv_nfold,
    # the leaves_only bootstrap re-estimates the random effects from the fold
    # models' out-of-fold residuals, as the point estimate does
    keep_fold_models = mse && bootstrap_refit == "leaves_only" && gbm_engine == "xgboost",
    em_maxit        = em_maxit,
    em_tol          = em_tol,
    ...
  )

  # ── 5b. Jensen-corrected back-transformation ─────────────────────────────────
  # E[g^{-1}(hat_t + e)] > g^{-1}(hat_t) for convex g^{-1}. We correct using
  # the MEGB-estimated error variance sigma2_e. Analytical formulae follow
  # expected_untransformed_mean() in point_estimation.R (exact
  # E[sin^2] formula for arcsin; exact lognormal formula for log/log.shift/poisson).
  sigma2_e <- tryCatch({
    err_sd <- megb_fit$megb_model$error_sd
    if (!is.null(err_sd) && !is.na(err_sd) && err_sd > 0) err_sd^2 else 0
  }, error = function(e) 0)

  # Jensen-corrected back-transformation (adds the random effect's conditional variance to sigma2_e)
  corrected_bt <- .megb_make_corrected_bt(sigma2_e, transformation,
    if (transformation == "log.shift") ls_lambda else NULL,
    arcsin_plain = arcsin_backtransform == "plain")
  # Marks the identity back-transform, the one case in which the leaves_only
  # bootstrap may aggregate from leaf values (povmap.leaves_only.leaf_aggregate).
  attr(corrected_bt, "identity") <- identical(transformation, "no")

  # ── 6. Back-transform unit predictions and aggregate to domain level ──────────
  unit_preds_pop      <- megb_fit$unit_preds_all       # columns: unit_preds, dom_name
  # Conditional variance of each domain's random effect (sigma_b^2 for unsampled domains,
  # the BLUP's prediction variance for sampled ones), from the same lmer fit that predicts
  # the random effect; added to sigma2_e in the back-transformation.
  re_cv_dom <- .megb_re_var(corrected_bt, megb_fit$megb_model$effect_model, all_levels)
  re_cv_row <- as.numeric(re_cv_dom[match(as.character(unit_preds_pop$dom_name), all_levels)])
  unit_preds_pop$hat  <- corrected_bt(unit_preds_pop$unit_preds, re_cv_row)

  if (!is.null(pop_weights)) {
    pw           <- fwk$pop_weights_vec                # same row order as pop_megb
    domain_means <- tapply(unit_preds_pop$hat * pw, unit_preds_pop$dom_name, sum) /
                    tapply(pw,                   unit_preds_pop$dom_name, sum)
  } else {
    domain_means <- tapply(unit_preds_pop$hat, unit_preds_pop$dom_name, mean)
  }

  ind <- data.frame(
    Domain = names(domain_means),
    Mean   = as.numeric(domain_means),
    stringsAsFactors = FALSE
  )

  # ── Benchmark key-alignment fix ───────────────────────────────────────────────
  # benchmark_xgb_level() positionally zips point_estim$ind$Mean against
  # unique(pop_data[[domains]]) (and the per-domain population weights). megb builds
  # `ind` via tapply(), which orders domains ALPHABETICALLY; pop_data is in
  # first-appearance order. When the two differ, every misplaced territory is raked
  # under the wrong benchmark-level target, silently corrupting the benchmarked
  # estimates (observed: ~90% of territories mis-assigned, with non-monotone
  # within-level redistribution and spurious zeros). xgb is unaffected because it
  # builds `ind` in pop_data order natively. Reorder `ind` to pop_data's domain
  # order so the positional contract holds; the assertion in benchmark_xgb_level
  # turns any residual mismatch into a loud error rather than silent corruption.
  .pop_dom_order <- as.character(unique(fwk$pop_data[[domains]]))
  ind <- ind[match(.pop_dom_order, as.character(ind$Domain)), , drop = FALSE]
  rownames(ind) <- NULL

  # predict_sampled = "crossfit": ind$Mean is cross-fitted. The same fit's
  # estimates with the full booster for every domain (Mean_full, what
  # predict_sampled = "full" returns) and with unsampled domains predicted by the
  # mean of the fold models (Mean_foldmean) are aggregated identically.
  crossfit_out <- NULL
  if (.fold_pred) {
    .agg <- function(u) {
      h <- corrected_bt(u, re_cv_row)
      m <- if (!is.null(pop_weights))
             tapply(h * fwk$pop_weights_vec, unit_preds_pop$dom_name, sum) /
             tapply(fwk$pop_weights_vec,     unit_preds_pop$dom_name, sum)
           else tapply(h, unit_preds_pop$dom_name, mean)
      as.numeric(m[match(ind$Domain, names(m))])
    }
    cf <- megb_fit$crossfit
    sampled <- ind$Domain %in% names(cf$fold_of_domain)
    crossfit_out <- list(
      ind = data.frame(Domain = ind$Domain, sampled = sampled,
                       fold = unname(cf$fold_of_domain[ind$Domain]),
                       Mean = ind$Mean, Mean_full = .agg(cf$unit_preds_full),
                       Mean_crossfit = .agg(cf$unit_preds_crossfit),
                       Mean_foldmean = .agg(cf$unit_preds_foldmean),
                       stringsAsFactors = FALSE),
      oof_max_abs_diff = cf$oof_max_abs_diff
    )
  }

  # ── 6b. Benchmarking of point estimates ───────────────────────────────────────
  ind_bench <- NULL
  if (!is.null(benchmark)) {
    bench_vals <- add_benchmark_megb(
      x               = ind$Mean,
      domain_labels   = ind$Domain,
      benchmark_level = benchmark_level,
      fwk             = fwk,
      fixed           = fixed,
      benchmark       = benchmark,
      benchmark_type  = benchmark_type
    )
    # Guard against silent benchmark failures (unsupported benchmark_type,
    # missing benchmark_level column, etc.) which previously surfaced as
    # all-empty Mean_bench / Var_bench / Lower_bench columns in write.excel.
    if (is.null(bench_vals) || length(bench_vals) != length(ind$Mean) ||
        all(is.na(bench_vals))) {
      stop(sprintf(
        "Benchmarking returned %s. Common causes: benchmark_type '%s' is not supported by the level-benchmarker (supported: 'ratio', 'ratio_bound', 'ratio_complement', 'logit_raking'); benchmark_level '%s' missing from smp_data or pop_data; or fixed[2] '%s' missing from smp_data.",
        if (is.null(bench_vals)) "NULL"
          else if (all(is.na(bench_vals))) "all-NA"
          else paste0("length ", length(bench_vals),
                      " (expected ", length(ind$Mean), ")"),
        as.character(benchmark_type),
        as.character(benchmark_level),
        as.character(fixed[[2]])
      ))
    }
    ind_bench <- data.frame(
      Domain = ind$Domain,
      Mean   = bench_vals,
      stringsAsFactors = FALSE
    )
  }

  # ── 6c. Build benchmark_fn closure for the bootstrap ─────────────────────────
  # When bench_target = "fixed" AND benchmark is character (internal), we
  # precompute the target ONCE from the observed survey using the same
  # survey-weighted-mean formula benchmark_xgb_level applies for internal
  # benchmarking. The closure then uses this fixed value instead of
  # recomputing from boot_smp_data$y_star each iteration, yielding a
  # conditional bench variance Var(theta_hat_bench | target = observed)
  # rather than the unconditional one that double-counts survey sampling
  # uncertainty on the target side.
  fixed_bm_val <- NULL
  if (!is.null(benchmark) && mse && bench_target == "fixed" &&
      !is.numeric(benchmark) && !is.data.frame(benchmark)) {
    bench_w_var <- if (!is.null(fwk$benchmark_weights)) fwk$benchmark_weights
                   else fwk$smp_weights
    if (!is.null(benchmark_level)) {
      grp_obs   <- fwk$smp_data[[benchmark_level]]
      y_obs     <- fwk$Y_smp
      w_obs     <- if (!is.null(bench_w_var)) fwk$smp_data[[bench_w_var]]
                   else rep(1, length(y_obs))
      lev_means <- tapply(seq_along(y_obs), grp_obs, function(idx)
        sum(y_obs[idx] * w_obs[idx]) / sum(w_obs[idx]))
      df <- data.frame(as.character(names(lev_means)),
                       as.numeric(lev_means),
                       stringsAsFactors = FALSE)
      names(df) <- c(benchmark_level, "Mean")
      fixed_bm_val <- df
    } else {
      # National-level internal benchmarking: scalar target.
      w_obs <- if (!is.null(bench_w_var)) fwk$smp_data[[bench_w_var]]
               else rep(1, length(fwk$Y_smp))
      fixed_bm_val <- c(Mean = sum(fwk$Y_smp * w_obs) / sum(w_obs))
    }
  }

  if (!is.null(benchmark) && mse) {
    # Use a lean captured environment: explicit values only, parent = package ns.
    # This prevents serialising the full megb() scope to each parallel worker.
    # domain_labels: the domain of each element of est_orig, in the order of est_orig. Required: the benchmark matches the
    # replicate's estimates to benchmark groups by these names (povmap <= 1.0.6 assumed the first-appearance order of
    # pop_data, which differs from the order of est_orig, so replicates were benchmarked under another domain's group).
    benchmark_fn <- function(est_orig, boot_smp_data, domain_labels) {
      if (missing(domain_labels) || is.null(domain_labels) || length(domain_labels) != length(est_orig))
        stop("benchmark_fn: domain_labels (one per element of est_orig) are required.")
      # Branch order: explicit user-supplied target → precomputed fixed target
      # (internal + bench_target="fixed") → recompute from bootstrap survey
      # (random target).
      # Random recompute MUST use the same survey-weighted mean as the point
      # benchmark (benchmark_xgb_level uses w = framework$smp_weights). An
      # unweighted mean would produce a different statistic and jitter far
      # more than the weighted target under skewed survey weights, widening
      # bench CIs spuriously.
      bm_val <- if (!is.null(attr(boot_smp_data, "bm_target"))) {
        attr(boot_smp_data, "bm_target")   # bench_target = "truth": this replicate's own target
      } else if (is.numeric(.bm) || is.data.frame(.bm)) {
        .bm
      } else if (!is.null(.fixed_bm)) {
        .fixed_bm
      } else if (is.null(.bm_level)) {
        w_obs <- if (!is.null(.bm_w_var)) boot_smp_data[[.bm_w_var]]
                 else rep(1, nrow(boot_smp_data))
        y_bt  <- .corr_bt(boot_smp_data$y_star)
        c(Mean = sum(y_bt * w_obs) / sum(w_obs))
      } else {
        y_bt  <- .corr_bt(boot_smp_data$y_star)
        grp   <- boot_smp_data[[.bm_level]]
        w_obs <- if (!is.null(.bm_w_var)) boot_smp_data[[.bm_w_var]]
                 else rep(1, length(y_bt))
        lev_means <- tapply(seq_along(y_bt), grp, function(idx)
          sum(y_bt[idx] * w_obs[idx]) / sum(w_obs[idx]))
        df <- data.frame(as.character(names(lev_means)),
                         as.numeric(lev_means),
                         stringsAsFactors = FALSE)
        names(df) <- c(.bm_level, "Mean")
        df
      }
      # benchmark_xgb_level discards non-numeric benchmark inputs and instead
      # recomputes from framework$smp_data (because of an `if (!is.numeric)`
      # branch upstream). For a data.frame target like our precomputed
      # .fixed_bm or the per-iteration random target above, that means our
      # value would be silently ignored. Convert to a named numeric vector,
      # which goes through the branch that actually uses the input.
      if (is.data.frame(bm_val) && !is.null(.bm_level)) {
        bm_val <- stats::setNames(as.numeric(bm_val[["Mean"]]),
                                  as.character(bm_val[[.bm_level]]))
      }
      .add_bm(
        x               = est_orig,
        domain_labels   = as.character(domain_labels),
        benchmark_level = .bm_level,
        fwk             = .fwk,
        fixed           = .fixed,
        benchmark       = bm_val,
        benchmark_type  = .bm_type
      )
    }
    # What the workers receive: the framework without the census columns benchmark_fn never reads, and the formula
    # without the environment it was written in (the whole megb() frame). Neither changes any number.
    fixed_lean <- fixed
    if (inherits(fixed_lean, "formula")) environment(fixed_lean) <- globalenv()
    environment(benchmark_fn) <- list2env(
      list(
        .bm       = benchmark,
        .bm_level = benchmark_level,
        .bm_type  = benchmark_type,
        .bm_w_var = if (!is.null(fwk$benchmark_weights)) fwk$benchmark_weights
                    else fwk$smp_weights,
        .fwk      = .megb_lean_framework(fwk, benchmark_level),
        .fixed    = fixed_lean,
        # y_star = gb + u* + e* already contains the unit error draw e*, like the observed y the
        # fixed target uses, so it takes the plain back-transform (sigma2_e = 0): no Jensen correction.
        .corr_bt  = .megb_make_corrected_bt(0, transformation,
                                            if (transformation == "log.shift") ls_lambda else NULL,
                                            arcsin_plain = arcsin_backtransform == "plain"),
        .add_bm   = add_benchmark_megb,
        .fixed_bm = fixed_bm_val
      ),
      parent = getNamespace("povmap")
    )
  } else {
    benchmark_fn <- NULL
  }

  # ── 6d. Parametric bootstrap MSE (now that corrected_bt/benchmark_fn are ready)
  # pop_data_proc / smp_data are one-hot encoded and lack the benchmark_level
  # column. Re-attach it from fwk$pop_data and fwk$smp_data so the bootstrap
  # samples (now derived from smp_data, not a pop subsample) carry it through
  # to benchmark_fn, which calls tapply(y_bt, boot_smp_data[[bm_level]], ...).
  boot_pop_data <- megb_fit$pop_data_proc
  boot_smp_data <- megb_fit$inp_smp_data$smp_data
  if (!is.null(benchmark_level) && mse) {
    if (!benchmark_level %in% colnames(boot_pop_data)) {
      bm_pop <- fwk$pop_data[[benchmark_level]]
      if (!is.null(bm_pop) && length(bm_pop) == nrow(boot_pop_data))
        boot_pop_data[[benchmark_level]] <- bm_pop
      else
        warning("benchmark_level '", benchmark_level,
                "' not found in pop_data or row count mismatch — ",
                "benchmarked bootstrap MSE will be skipped.")
    }
    if (!benchmark_level %in% colnames(boot_smp_data)) {
      bm_smp <- fwk$smp_data[[benchmark_level]]
      if (!is.null(bm_smp) && length(bm_smp) == nrow(boot_smp_data))
        boot_smp_data[[benchmark_level]] <- bm_smp
      else
        warning("benchmark_level '", benchmark_level,
                "' not found in smp_data or row count mismatch — ",
                "benchmarked bootstrap MSE will be skipped.")
    }
  }

  # bench_target = "truth": the group of each domain, its population total (the weight benchmark_xgb_level aggregates the
  # estimate with) and the B x G draws of the target error, made once before the bootstrap from a stream of their own.
  bench_star <- NULL
  if (bench_target == "truth" && mse) {
    dom_lab   <- as.character(fwk$pop_data[[domains]])
    dom_u     <- sort(unique(dom_lab))   # sorted, not first-appearance: the draws are assigned to groups in this order, so they must not depend on the row order of pop_data (as in .megb_perturb_bench)
    dom_group <- stats::setNames(as.character(fwk$pop_data[[benchmark_level]][match(dom_u, dom_lab)]), dom_u)
    dom_w     <- tapply(as.numeric(fwk$pop_data[[fwk$pop_weights]]), dom_lab, sum)
    bm_draws  <- .megb_bench_draws(benchmark_target_se, unique(dom_group[!is.na(dom_group)]), B, seed,
                                   allow_zero_benchmark_se, "bench_target = \"truth\"")$draws
    bm_bounds <- .megb_bench_bounds(transformation)
    bench_star <- list(dom_group = dom_group, dom_w = stats::setNames(as.numeric(dom_w), names(dom_w)),
                       draws = bm_draws, lo = bm_bounds[["lo"]], hi = bm_bounds[["hi"]])
  }

  mse_estimated <- NULL
  if (mse) {
    message("Bootstrap with ", B, " iterations has started")
    # Attach weights as a column on boot_smp_data to guarantee row alignment.
    # fwk$smp_weights_vec was built from fwk$smp_data, but boot_smp_data
    # (= megb_fit$inp_smp_data$smp_data) may have a different row count or
    # ordering after framework + megb_em processing. Pulling weights from a
    # column ensures they always match the smp_data rows downstream consumers
    # iterate over.
    if (!is.null(fwk$smp_weights_vec)) {
      if (length(fwk$smp_weights_vec) == nrow(boot_smp_data)) {
        boot_smp_data$.megb_w <- as.numeric(fwk$smp_weights_vec)
      } else {
        warning("smp_weights_vec length (", length(fwk$smp_weights_vec),
                ") does not match boot_smp_data row count (",
                nrow(boot_smp_data), "). Falling back to unweighted bootstrap.")
        boot_smp_data$.megb_w <- NULL
      }
    }
    mse_estimated <- mse_megb(
      Y                      = megb_fit$inp_smp_data$target_var,
      X                      = megb_fit$X_proc,
      dom_name               = domains,
      smp_data               = boot_smp_data,
      model                  = megb_fit$megb_model,
      error_sd               = megb_fit$megb_model$error_sd,
      pop_data               = boot_pop_data,
      B                      = B,
      initial_random_effects = 0,
      ErrorTolerance         = em_tol,
      MaxIterations          = em_maxit,
      cov_names              = megb_fit$cov_names_proc,
      gradient_params        = gradient_params,
      formula_random_effects = megb_fit$formula_random_effects,
      bootstrap_cores        = bootstrap_cores,
      seed                   = seed,
      gbm_engine             = gbm_engine,
      unit_pred_smp          = megb_fit$unit_pred_smp,
      gb_smp                 = megb_fit$gb_smp,
      unit_preds             = megb_fit$unit_preds_all,
      bootstrap_refit        = bootstrap_refit,
      corrected_bt           = corrected_bt,
      benchmark_fn           = benchmark_fn,
      smp_weights_col        = ".megb_w",
      weightedBS             = weightedBS,
      # Same per-cell population weights the point estimate aggregates with
      # (megb.R: domain_means = sum(hat*pw)/sum(pw) over pop cells). Passed so the
      # bootstrap's per-cell-then-aggregate domain means use identical weighting;
      # NULL reproduces the historical unweighted bootstrap aggregation.
      pop_weights_vec        = if (!is.null(pop_weights)) fwk$pop_weights_vec else NULL,
      cv_folds               = cv_folds,
      early_stopping         = early_stopping,
      cv_nfold               = cv_nfold,
      fold_fit               = .megb_fold_fit(megb_fit$megb_model),
      predict_sampled        = predict_sampled,
      predict_unsampled      = predict_unsampled,
      bench_star             = bench_star,
      boot_area_var          = boot_area_var,
      boot_synthetic         = boot_synthetic
    )
  }

  # ── 6e. Benchmark-target perturbation (perturb_benchmark = TRUE) ──────────────
  # As in xgb: after the bootstrap, add the pre-drawn target noise to the benchmarked replicates on the
  # model's transformed scale; the truth side is not perturbed.
  bm_perturbation <- NULL
  if (perturb_benchmark && !is.null(mse_estimated$tau_b_bench)) {
    mse_estimated <- .megb_perturb_bench(
      mse_estimated, benchmark_target_se = benchmark_target_se,
      allow_zero = allow_zero_benchmark_se, B = B, seed = seed,
      dom_group = fwk$pop_data[[benchmark_level]][match(as.character(mse_estimated$MSE_estimates[[1]]),
                                                         as.character(fwk$pop_data[[domains]]))],
      transformation = transformation, transform_outcome = transform_outcome,
      back_transform_outcome = back_transform_outcome)
    bm_perturbation <- mse_estimated$benchmark_perturbation
  }

  # ── 7. Variance / MSE and CI ─────────────────────────────────────────────────
  # mse_type = "var" (default, mirrors xgb): empirical bootstrap variance and
  #   quantile-based CIs in original space. Benchmarked CIs use the centred
  #   residual (boot_bench - boot_truth) recentred at the point benchmark, as
  #   in xgb.R.
  # mse_type = "mse": prediction MSE in original space. Unbenchmarked MSE is
  #   computed in transformed space from the bootstrap and back-mapped via the
  #   delta method; benchmarked MSE is computed directly in original space.
  #   CIs use a normal approximation around the point estimate.
  var_df       <- NULL
  ci_df        <- NULL
  var_bench_df <- NULL
  ci_bench_df  <- NULL

  if (mse && !is.null(mse_estimated)) {

    alpha <- 1 - conf_level
    z_val <- qnorm(1 - alpha / 2)

    # Helper: clamp CIs to the valid outcome support for bounded transformations.
    clamp_ci <- function(df) {
      if (transformation %in% c("arcsin", "logistic")) {
        df$Lower <- pmax(df$Lower, 0); df$Upper <- pmin(df$Upper, 1)
      } else if (transformation %in% c("poisson", "log", "log.shift")) {
        df$Lower <- pmax(df$Lower, 0)
      }
      df
    }

    boot_domains <- as.character(mse_estimated$domains)

    if (mse_type == "var") {

      # ── 7a-var. Unbenchmarked: centred residuals of (boot - boot_truth) ──
      # For a mixed-effects estimator the target of inference is the realised
      # domain mean (which depends on u_d), so the natural uncertainty summary
      # is the prediction error against the bootstrap truth, not the variance
      # of the bootstrap estimates alone. This is symmetric with how the
      # benchmarked CI is computed below and with the "mse" path.
      tau_b_orig <- mse_estimated$tau_b_orig         # D × B' (original space)
      truth_u    <- mse_estimated$tau_star_orig_unb  # truth restricted to same B'

      resid_u <- tau_b_orig - truth_u
      EB <- .megb_error_interval(resid_u, alpha)

      ind$Domain <- as.character(ind$Domain)
      idx        <- match(ind$Domain, boot_domains)

      var_df <- data.frame(Domain = ind$Domain, Mean = EB$var[idx],
                           stringsAsFactors = FALSE)

      ci_df <- data.frame(
        Domain = ind$Domain,
        Lower  = ind$Mean - EB$q_hi[idx],
        Upper  = ind$Mean - EB$q_lo[idx],
        stringsAsFactors = FALSE
      )
      ci_df <- clamp_ci(ci_df)

      # ── 7b-var. Benchmarked: var/quantiles of (boot_bench - boot_truth) ──
      # Both matrices are restricted to the same kept iterations by mse_megb,
      # so column counts match without further index gymnastics.
      if (!is.null(benchmark) && !is.null(mse_estimated$tau_b_bench) &&
          !is.null(mse_estimated$tau_star_orig_bench)) {
        tau_b_bench <- mse_estimated$tau_b_bench
        truth_bench <- mse_estimated$tau_star_orig_bench
        if (ncol(tau_b_bench) > 0L) {
          EBb <- .megb_error_interval(tau_b_bench - truth_bench, alpha)
          var_bench_vec <- EBb$var

          ind_bench$Domain <- as.character(ind_bench$Domain)
          idx_b <- match(ind_bench$Domain, boot_domains)

          var_bench_df <- data.frame(Domain = ind_bench$Domain,
                                     Mean   = var_bench_vec[idx_b],
                                     stringsAsFactors = FALSE)
          ci_bench_df <- data.frame(
            Domain = ind_bench$Domain,
            Lower  = ind_bench$Mean - EBb$q_hi[idx_b],
            Upper  = ind_bench$Mean - EBb$q_lo[idx_b],
            stringsAsFactors = FALSE
          )
          ci_bench_df <- clamp_ci(ci_bench_df)
        }
      }

    } else {  # mse_type == "mse"

      # ── 7a-mse. Unbenchmarked prediction MSE (delta-method back-map) ──
      mse_raw              <- mse_estimated$MSE_estimates
      colnames(mse_raw)[1] <- "Domain"
      mse_raw$Domain       <- as.character(mse_raw$Domain)
      colnames(mse_raw)[2] <- "MSE_t"

      mean_t        <- as.data.frame(megb_fit$Indicators)
      colnames(mean_t)[colnames(mean_t) == "dom_name"] <- "Domain"
      colnames(mean_t)[colnames(mean_t) == "Mean"]     <- "Mean_t"
      mean_t$Domain <- as.character(mean_t$Domain)

      ind$Domain <- as.character(ind$Domain)
      merged     <- merge(ind,    mse_raw, by = "Domain")
      merged     <- merge(merged, mean_t,  by = "Domain")

      delta_factor <- switch(transformation,
        "no"        = rep(1, nrow(merged)),
        "log"       = exp(merged$Mean_t)^2,
        "log.shift" = exp(merged$Mean_t)^2,
        "arcsin"    = sin(2 * merged$Mean_t)^2,
        "logistic"  = { p <- back_transform_outcome(merged$Mean_t); (p * (1 - p))^2 },
        "poisson"   = exp(merged$Mean_t)^2
      )
      merged$MSE_orig <- merged$MSE_t * delta_factor

      var_df <- data.frame(Domain = merged$Domain, Mean = merged$MSE_orig,
                           stringsAsFactors = FALSE)
      ci_df <- data.frame(
        Domain = merged$Domain,
        Lower  = merged$Mean - z_val * sqrt(pmax(merged$MSE_orig, 0)),
        Upper  = merged$Mean + z_val * sqrt(pmax(merged$MSE_orig, 0)),
        stringsAsFactors = FALSE
      )
      ci_df <- clamp_ci(ci_df)
      ind   <- data.frame(Domain = merged$Domain, Mean = merged$Mean,
                          stringsAsFactors = FALSE)

      # ── 7b-mse. Benchmarked prediction MSE (already in original space) ──
      if (!is.null(benchmark) && !is.null(mse_estimated$MSE_bench_estimates)) {
        mse_bench_raw              <- mse_estimated$MSE_bench_estimates
        colnames(mse_bench_raw)[1] <- "Domain"
        mse_bench_raw$Domain       <- as.character(mse_bench_raw$Domain)
        colnames(mse_bench_raw)[2] <- "MSE_bench"

        ind_bench$Domain <- as.character(ind_bench$Domain)
        merged_b         <- merge(ind_bench, mse_bench_raw, by = "Domain")

        var_bench_df <- data.frame(Domain = merged_b$Domain, Mean = merged_b$MSE_bench,
                                   stringsAsFactors = FALSE)
        ci_bench_df <- data.frame(
          Domain = merged_b$Domain,
          Lower  = merged_b$Mean - z_val * sqrt(pmax(merged_b$MSE_bench, 0)),
          Upper  = merged_b$Mean + z_val * sqrt(pmax(merged_b$MSE_bench, 0)),
          stringsAsFactors = FALSE
        )
        ci_bench_df <- clamp_ci(ci_bench_df)
        ind_bench   <- data.frame(Domain = merged_b$Domain, Mean = merged_b$Mean,
                                  stringsAsFactors = FALSE)
      }
    }
  }

  # Unperturbed benchmarked variance (var) or MSE (mse), kept next to the perturbed one as in xgb.
  var_bench_unpert_df <- NULL
  if (!is.null(mse_estimated$tau_b_bench_unperturbed)) {
    err0 <- mse_estimated$tau_b_bench_unperturbed - mse_estimated$tau_star_orig_bench
    v0   <- if (mse_type == "var") .megb_error_interval(err0, alpha)$var else rowMeans(err0^2, na.rm = TRUE)
    var_bench_unpert_df <- data.frame(Domain = as.character(mse_estimated$MSE_estimates[[1]]),
                                      Var_bench_unperturbed = as.numeric(v0), stringsAsFactors = FALSE)
  }

  # ── 8. Sample-level predictions (back-transformed) ───────────────────────────
  unit_pred_smp_bt <- corrected_bt(megb_fit$unit_pred_smp)
  yhat <- data.frame(
    obs_id = seq_along(unit_pred_smp_bt),
    hat    = as.numeric(unit_pred_smp_bt)
  )

  # Store in framework for summary use
  fwk$unit_pred_smp <- as.numeric(unit_pred_smp_bt)

  # ── 9. Model reference ────────────────────────────────────────────────────────
  # For xgboost: expose the booster so add_model_megb can compute importance
  if (gbm_engine == "xgboost") {
    model_ref <- megb_fit$megb_model$boosting
  } else {
    model_ref <- megb_fit$megb_model
  }

  # ── 9b. Align result structure with xgb so shared write.excel/estimators logic works ──
  # xgb stores Mean_bench as a column in $ind, Var_bench in $var, Lower/Upper_bench in $CI.
  if (!is.null(ind_bench)) {
    ind <- merge(
      ind,
      data.frame(Domain = ind_bench$Domain, Mean_bench = ind_bench$Mean,
                 stringsAsFactors = FALSE),
      by = "Domain", all.x = TRUE
    )
  }
  if (!is.null(var_df) && !is.null(var_bench_df)) {
    var_df <- merge(
      var_df,
      data.frame(Domain = var_bench_df$Domain, Var_bench = var_bench_df$Mean,
                 stringsAsFactors = FALSE),
      by = "Domain", all.x = TRUE
    )
  }
  if (!is.null(ci_df) && !is.null(ci_bench_df)) {
    ci_bench_renamed <- ci_bench_df
    names(ci_bench_renamed)[names(ci_bench_renamed) == "Lower"] <- "Lower_bench"
    names(ci_bench_renamed)[names(ci_bench_renamed) == "Upper"] <- "Upper_bench"
    ci_df <- merge(ci_df, ci_bench_renamed, by = "Domain", all.x = TRUE)
  }

  # ── 9c. Bootstrap diagnostics ────────────────────────────────────────────────
  # Surface boot_ran_eff_sd / boot_error_sd alongside their non-bootstrap
  # counterparts so the caller can diagnose pathologies like "the EM in each
  # bootstrap iteration collapses to ran_eff_sd ~ 0", which would inflate the
  # prediction-MSE summary by removing the refit's ability to recover u_d*.
  boot_diag <- NULL
  if (!is.null(mse_estimated)) {
    qsum <- function(x) {
      x <- x[is.finite(x)]
      if (length(x) == 0L) return(rep(NA_real_, 5L))
      as.numeric(quantile(x, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE))
    }
    boot_diag <- list(
      ran_eff_sd_orig    = megb_fit$megb_model$ran_eff_sd,
      error_sd_orig      = megb_fit$megb_model$error_sd,
      boot_ran_eff_sd    = mse_estimated$boot_ran_eff_sd_boot,
      boot_error_sd      = mse_estimated$boot_error_sd,
      boot_em_iterations = mse_estimated$boot_em_iterations,
      boot_ran_eff_sd_q  = qsum(mse_estimated$boot_ran_eff_sd_boot),
      boot_error_sd_q    = qsum(mse_estimated$boot_error_sd)
    )
    # boot_area_var != "legacy": mode, the moment estimator sigma_b2(MM) and the pool sizes (absent for "legacy")
    if (!identical(boot_area_var, "legacy")) boot_diag$boot_area_var <- mse_estimated$boot_area_mm
    # boot_synthetic = TRUE: the point fit's synthetic area prediction, the replicates' (reporting scale, unbenchmarked), their truth and the MSE of the
    # synthetic prediction per area; the input of megb_level_calibration() (MEGB_BOOT_FIX_PREREG_v8_final.md, option A)
    if (boot_synthetic && !is.null(mse_estimated$synth_boot_orig)) {
      sy_dom <- as.character(mse_estimated$domains)
      re_syn <- .megb_re_var(corrected_bt, megb_fit$megb_model$effect_model, rep(".synthetic.", nrow(unit_preds_pop)))
      h_syn  <- corrected_bt(unit_preds_pop$gb_pop, re_syn)
      sy_pt  <- if (!is.null(pop_weights)) tapply(h_syn * pw, unit_preds_pop$dom_name, sum) / tapply(pw, unit_preds_pop$dom_name, sum)
                else tapply(h_syn, unit_preds_pop$dom_name, mean)
      sy_truth <- mse_estimated$tau_star_orig_unb
      boot_diag$synthetic <- list(
        domains   = sy_dom,
        sampled   = sy_dom %in% as.character(unique(fwk$smp_data[[domains]])),
        point     = as.numeric(sy_pt[match(sy_dom, names(sy_pt))]),
        boot      = mse_estimated$synth_boot_orig,
        truth     = sy_truth,
        mse_synth = rowMeans((mse_estimated$synth_boot_orig - sy_truth)^2, na.rm = TRUE))
    }
    message(sprintf(
      "Bootstrap diagnostics: original ran_eff_sd = %.4f, error_sd = %.4f.\n  Bootstrap ran_eff_sd quartiles (min,Q1,Q2,Q3,max) = %s\n  Bootstrap error_sd   quartiles                  = %s\n  Frac bootstrap iters with ran_eff_sd < 1e-4 = %.2f",
      megb_fit$megb_model$ran_eff_sd %||% NA_real_,
      megb_fit$megb_model$error_sd   %||% NA_real_,
      paste(sprintf("%.4f", boot_diag$boot_ran_eff_sd_q), collapse = ", "),
      paste(sprintf("%.4f", boot_diag$boot_error_sd_q),   collapse = ", "),
      mean(boot_diag$boot_ran_eff_sd < 1e-4, na.rm = TRUE)
    ))
  }

  # ── 10. Assemble result ───────────────────────────────────────────────────────
  result <- list(
    ind            = ind,
    var            = var_df,
    CI             = ci_df,
    ind_bench      = ind_bench,
    var_bench      = var_bench_df,
    CI_bench       = ci_bench_df,
    yhat           = yhat,
    model          = model_ref,
    megb_model     = megb_fit$megb_model,
    gbm_engine     = gbm_engine,
    smp_data       = smp_data,
    out_call       = out_call,
    transformation = transformation,
    framework      = fwk,
    boot_diag      = boot_diag,
    predict_sampled = predict_sampled,
    predict_unsampled = predict_unsampled,
    cv_nfold       = cv_nfold,
    crossfit       = crossfit_out
  )

  if (perturb_benchmark) {
    result$Var_bench_unperturbed <- var_bench_unpert_df
    result$benchmark_perturbation <- bm_perturbation
  }
  if (identical(transformation, "arcsin") && arcsin_backtransform == "plain") result$arcsin_backtransform <- "plain"
  if (!is.null(bench_star)) {
    result$bench_target <- "truth"
    result$benchmark_target_draws <- bench_star$draws   # B x G draws of the target error (N(0, SE_g^2)), by benchmark group
  }

  class(result) <- c("megb", "xgb", "povmap")
  return(result)
}
