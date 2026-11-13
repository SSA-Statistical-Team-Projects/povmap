# DRC follow-up retune: edge picks of the CV3 run, same protocol, wider bounds only where a pick sat on an edge.
# Extends drc_retune_repeatedcv.R. NOT TO BE RUN while another job holds the machine (see STAGE2_PLAN.md for the run order and worker plan).
#
# Why these indicators (CV3 picks on a bound, drc_retune_repeatedcv.csv):
#   petitrum      lambda 9055 against an upper bound of 10000          -> lambda upper bound 1e5 (log scale)
#   pcexp, lfp, telpor   max_depth 10 against an upper bound of 10    -> max_depth upper bound 15
#   school_12_17, porc   max_depth 2 against a LOWER bound of 2       -> max_depth lower bound 1 (stumps)
# Other picks at an edge of the CV3 bounds (reported in the output csv for all six, and listed in STAGE2_PLAN.md for the other ten;
# NOT widened here because the instruction names only the cases above): colsample_bytree at its 0.3 lower edge (lfp, telpor, school_12_17,
# porc), min_child_weight 1.09-1.23 against a lower bound of 1 (nivie, theil, access_water), nround 1492 of 1500 (access_elec), eta 0.266 of 0.3
# (volail), subsample 0.522 of 0.5 (tv).
#
# Candidates. The CV3 candidate table (150 draws, seed 20261010, from repeatedcv/detail/<indicator>.rds) is reused, so nothing is scored on
# seed 20261010 to obtain draws:
#   * lambda widening: the sampler draws log-uniformly from the same uniform stream whatever the bound, so the widened draws are an exact
#     rescaling of the old ones, lambda_new = lo * (lambda_old / lo)^(log(hi_new / lo) / log(hi_old / lo)); every other parameter is
#     identical (checked against the package sampler, .xgb_sample_bounds, in the log).
#   * max_depth widening: depth is drawn with sample(lo:hi, n, replace = TRUE) as the SECOND parameter of the stream, and sample() consumes a
#     different number of random numbers when the range changes, so re-running the sampler with a wider depth range changes EVERY later
#     parameter's draws (checked: eta is the only unchanged column). To keep all other parameters as they were, the depth column alone is
#     redrawn here, from its own stream (set.seed(SEED + 15), sample(lo:hi, 150, TRUE)), integer-uniform on the widened range; every other
#     column is copied from the CV3 table. The candidate set therefore differs from CV3's only in depth (and, for petitrum, lambda).
# Scoring and selection are exactly the CV3 protocol: fold seeds 123, 124, 125 (and 20261010), response scale, fit seed 123, nthread 1,
# paired one-SE set (floor 1, no cap), least-overfit on the area (domain) pool, repeated-CV rule xgb_top_configs(list of the three tunes),
# area ratio computed on each seed's folds and averaged. The CV3 pick and set B are scored alongside (one-row tunes on the same folds).
# COMBINED MODE (ALL=1): the rescore of all 16 indicators on SCORE_SCALE and the follow-up in one run. The six widened indicators use their widened
# candidate tables, the other ten the CV3 tables unchanged, so each candidate is fitted once. Only fold seeds 123, 124, 125 are scored (WITH_DRAW_SEED=1
# adds 20261010). For the ten unchanged indicators the overfit ratios already computed in CV3 (they do not depend on the score scale) are reused for the
# rows that have them and only the missing rows are fitted.
# Environment: INDICATORS (default the six above), ALL (1 = all 16), NW (workers), OUT, SEED (20261010), TAG, FINAL,
#   SCORE_SCALE ("transformed" (default), "response" or "training"; the CV3 run used "response"), WITH_DRAW_SEED (0/1).
suppressMessages(library(povmap))
DATA <- Sys.getenv("DRC_DATA", "H:/david/SAE/DRC/data")
BASE <- file.path(DATA, "candidate_estimates", "drc_retune_20261010")
OUT <- Sys.getenv("OUT", file.path(BASE, "followup"))
NW <- as.integer(Sys.getenv("NW", "6")); SEED <- as.integer(Sys.getenv("SEED", "20261010")); FIT_SEED <- 123L
SEEDS <- c(123L, 124L, 125L)
ALLSEEDS <- if (Sys.getenv("WITH_DRAW_SEED", "0") == "1") c(SEED, SEEDS) else SEEDS
SCORE_SCALE <- Sys.getenv("SCORE_SCALE", "transformed"); stopifnot(SCORE_SCALE %in% c("transformed", "response", "training"))
dir.create(file.path(OUT, "detail"), showWarnings = FALSE, recursive = TRUE)
survey <- readRDS(file.path(DATA, "survey_with_features")); census <- readRDS(file.path(DATA, "census_xgb_agepop"))
X <- colnames(survey)[c(which(colnames(survey) == "population"), which(colnames(survey) == "ntl"):which(colnames(survey) == "FATALITIES"),
                        which(colnames(survey) == "population_territory"), which(colnames(survey) == "ntl_territory"):length(colnames(survey)))]
X <- intersect(setdiff(X, c("popwt", "hhwt", "pop_weight", "wt_5_11", "wt_12_17", "wt_15_64")), names(census))
prov_map <- readRDS(file.path(DATA, "prov_fe_map.rds")); prov_cols <- names(prov_map)
for (pc in prov_cols) survey[[pc]] <- as.integer(as.character(survey$NAM_1) == prov_map[[pc]])
X <- c(X, prov_cols)
cfg <- list(nivie = c("arcsin", "popwt", NA), theil = c("log", "popwt", NA), pcexp = c("log", "popwt", NA), access_elec = c("arcsin", "popwt", NA), access_water = c("arcsin", "popwt", NA),
            lfp = c("arcsin", "wt_15_64", "lfp_var"), school_5_11 = c("arcsin", "wt_5_11", "school_5_11_var"), school_12_17 = c("arcsin", "wt_12_17", "school_12_17_var"),
            telpor = c("arcsin", "popwt", NA), grosrum = c("poisson", "hhwt", NA), petitrum = c("poisson", "hhwt", NA), porc = c("poisson", "hhwt", NA), lapin = c("poisson", "hhwt", NA),
            volail = c("poisson", "hhwt", NA), car = c("arcsin", "popwt", NA), tv = c("arcsin", "popwt", NA))
CV3 <- list(eta = c(0.01, 0.3), max_depth = c(2, 10), min_child_weight = c(1, 1000), subsample = c(0.5, 1), colsample_bytree = c(0.3, 1),
            lambda = c(0.01, 10000), nround = c(50, 1500))
WIDEN <- list(petitrum = list(lambda = c(0.01, 1e5)), pcexp = list(max_depth = c(2, 15)), lfp = list(max_depth = c(2, 15)), telpor = list(max_depth = c(2, 15)),
              school_12_17 = list(max_depth = c(1, 10)), porc = list(max_depth = c(1, 10)))
todo <- if (nzchar(Sys.getenv("INDICATORS"))) strsplit(Sys.getenv("INDICATORS"), ",")[[1]] else if (Sys.getenv("ALL") == "1") names(cfg) else names(WIDEN)
KEYS <- c("nround", "max_depth", "colsample_bytree", "colsample_bylevel", "colsample_bynode", "subsample", "min_child_weight", "eta", "gamma", "max_delta_step", "lambda", "alpha")
SHOW <- c("nround", "max_depth", "eta", "min_child_weight", "lambda", "subsample", "colsample_bytree")
cv3 <- read.csv(file.path(BASE, "drc_retune_repeatedcv.csv"), stringsAsFactors = FALSE)
setB <- readRDS(file.path(DATA, "candidate_estimates", "least_overfit_20261008", "xgb_all_tunes_least_overfit"))
retry <- function(expr, tries = 6) {
  e <- substitute(expr); pf <- parent.frame()
  for (i in seq_len(tries)) {
    r <- tryCatch(eval(e, pf), error = function(err) if (grepl("serverSocket|port", conditionMessage(err)) && i < tries) NULL else stop(err))
    if (!is.null(r)) return(r)
    Sys.sleep(5)
  }
}
rows_of <- function(tn, n) attr(xgb_top_configs(tn, rule = "one_se_ratio", ratio = rep(1, n), max_ratio = 2, cap = n, ratio_cap = NULL), "rows")
bound_flags <- function(v, bounds) {
  pos <- povmap:::.xgb_bounds_position(bounds, as.list(v)); f <- pos[pos$at_boundary, ]
  if (!nrow(f)) return("")
  paste(sprintf("%s=%s (%s bound %s)", f$parameter, signif(f$selected, 4), ifelse(f$position < 0.5, "lower", "upper"), ifelse(f$position < 0.5, f$lower, f$upper)), collapse = "; ")
}
res <- list()
for (ind in todo) {
  tr <- cfg[[ind]][1]; wc <- cfg[[ind]][2]; hc <- cfg[[ind]][3]; hv <- if (is.na(hc)) NULL else hc
  bounds <- CV3; widen <- if (ind %in% names(WIDEN)) WIDEN[[ind]] else list(); for (p in names(widen)) bounds[[p]] <- widen[[p]]
  old <- readRDS(file.path(BASE, "repeatedcv", "detail", paste0(ind, ".rds")))$candidates; n <- nrow(old)
  cand <- old[, KEYS]
  if ("lambda" %in% names(widen)) {
    lo <- CV3$lambda[1]; cand$lambda <- lo * (old$lambda / lo)^(log(bounds$lambda[2] / lo) / log(CV3$lambda[2] / lo))
    ## (checked in the planning step against .xgb_sample_bounds: the rescaling reproduces the sampler exactly, all other columns unchanged)
  }
  if ("max_depth" %in% names(widen)) {
    had <- exists(".Random.seed", envir = globalenv()); if (had) sv <- get(".Random.seed", envir = globalenv())
    set.seed(SEED + 15L); cand$max_depth <- sample(seq.int(bounds$max_depth[1], bounds$max_depth[2]), n, replace = TRUE)
    if (had) assign(".Random.seed", sv, envir = globalenv())
  }
  smp <- survey[!is.na(survey[[ind]]), unique(c(ind, X, wc, hv, "ADM2CD_c"))]
  fx <- stats::reformulate(X, response = ind)
  common <- list(fx, smp_data = smp, smp_weights = wc, domains = "ADM2CD_c", transformation = tr, folds = 10, fit_seed = FIT_SEED,
                 nthread = 1, variance_y = hv, score_scale = SCORE_SCALE, cpus = NW, verbose = FALSE)
  tn <- list(); for (s in ALLSEEDS) tn[[paste0("s", s)]] <- retry(suppressWarnings(do.call(xgb_tune, c(common, list(seed = s, grid = cand)))))
  rep_tn <- tn[paste0("s", SEEDS)]
  sets <- c(lapply(tn, rows_of, n = n), list(repeated = rows_of(unname(rep_tn), n)))
  U <- sort(unique(unlist(sets)))
  ## ratios do not depend on the score scale: reuse the CV3 ones for unchanged candidate tables, fit only the missing rows
  cv3_old <- if (length(widen) == 0L) readRDS(file.path(BASE, "repeatedcv", "detail", paste0(ind, ".rds")))$tunes else NULL
  n_ratio_fits <- 0L
  for (s in ALLSEEDS) {
    ov <- data.frame(row = seq_len(n), ratio_total = NA_real_, ratio_domain = NA_real_, ratio_within_domain = NA_real_)
    prev <- if (!is.null(cv3_old)) cv3_old[[paste0("s", s)]]$overfit else NULL
    if (!is.null(prev)) ov[, c("ratio_total", "ratio_domain", "ratio_within_domain")] <- prev[, c("ratio_total", "ratio_domain", "ratio_within_domain")]
    need <- U[is.na(ov$ratio_domain[U])]
    if (length(need)) {
      rt <- retry(xgb_overfit_ratio(fx, smp, cand, smp_weights = wc, domains = "ADM2CD_c", transformation = tr, folds = 10, variance_y = hv,
                                    rows = need, seed = s, fit_seed = FIT_SEED, cpus = NW))
      ov[rt$row, c("ratio_total", "ratio_domain", "ratio_within_domain")] <- rt[, c("ratio_total", "ratio_domain", "ratio_within_domain")]
      n_ratio_fits <- n_ratio_fits + nrow(rt)
    }
    tn[[paste0("s", s)]]$overfit <- ov
  }
  rep_tn <- tn[paste0("s", SEEDS)]
  pk <- function(t) attr(xgb_top_configs(t, rule = "least_overfit", ratio_pool = "domain", ratio_cap = NULL), "rows")
  picks <- c(lapply(tn, pk), list(repeated = pk(unname(rep_tn))))
  loss_rep <- rowMeans(vapply(rep_tn, function(t) t$candidates$mse, numeric(n)))
  ratio_rep <- list(area = rowMeans(vapply(rep_tn, function(t) t$overfit$ratio_domain, numeric(n))),
                    total = rowMeans(vapply(rep_tn, function(t) t$overfit$ratio_total, numeric(n))),
                    within = rowMeans(vapply(rep_tn, function(t) t$overfit$ratio_within_domain, numeric(n))))
  p <- picks$repeated; b <- which.min(loss_rep)
  ## the CV3 pick and set B scored on the same folds
  c3 <- cv3[cv3$indicator == ind, ]; g3 <- old[c3$pick_row_repeated, KEYS]
  sc <- function(g) { tt <- lapply(SEEDS, function(s) retry(suppressWarnings(do.call(xgb_tune, c(common[names(common) != "cpus"], list(cpus = NW, seed = s, grid = g, overfit_ratio = TRUE))))))
                      c(loss = mean(vapply(tt, function(t) t$candidates$mse[1], 1)), area = mean(vapply(tt, function(t) t$overfit$ratio_domain[1], 1)),
                        total = mean(vapply(tt, function(t) t$overfit$ratio_total[1], 1)), within = mean(vapply(tt, function(t) t$overfit$ratio_within_domain[1], 1))) }
  s3 <- sc(g3)
  r <- data.frame(indicator = ind, widening = if (length(widen)) paste(names(widen), vapply(widen, function(v) paste(v, collapse = "-"), ""), collapse = "; ") else "none", score_scale = SCORE_SCALE, ratio_rows_fitted = n_ratio_fits,
                  set_size_repeated = length(sets$repeated), pick_row = p, cand[p, SHOW],
                  pick_cv_loss = loss_rep[p], best_cv_loss = loss_rep[b], pick_ratio_area = ratio_rep$area[p], pick_ratio_total = ratio_rep$total[p], pick_ratio_within = ratio_rep$within[p],
                  pick_row_s20261010 = if (is.null(picks$s20261010)) NA_integer_ else picks$s20261010, pick_row_s123 = picks$s123, pick_row_s124 = picks$s124, pick_row_s125 = picks$s125,
                  cv3_pick_row = c3$pick_row_repeated, cv3_nround = c3$nround, cv3_max_depth = c3$max_depth, cv3_eta = c3$eta, cv3_min_child_weight = c3$min_child_weight,
                  cv3_lambda = c3$lambda, cv3_subsample = c3$subsample, cv3_colsample_bytree = c3$colsample_bytree,
                  cv3_pick_cv_loss_on_this_protocol = s3[["loss"]], cv3_pick_ratio_area_on_this_protocol = s3[["area"]],
                  new_pick_equals_cv3_pick = isTRUE(all(abs(as.numeric(cand[p, SHOW]) - as.numeric(old[c3$pick_row_repeated, SHOW])) <= 1e-9 * pmax(1, abs(as.numeric(old[c3$pick_row_repeated, SHOW]))))),
                  setB_nround = c3$setB_nround, setB_max_depth = c3$setB_max_depth, setB_eta = c3$setB_eta, setB_min_child_weight = c3$setB_min_child_weight,
                  setB_lambda = c3$setB_lambda, setB_subsample = c3$setB_subsample, setB_colsample_bytree = c3$setB_colsample_bytree,
                  setB_cv_loss = c3$setB_cv_loss_repeated, setB_ratio_area = c3$setB_ratio_area,
                  boundary_flags_new_pick = bound_flags(cand[p, ], bounds), boundary_flags_new_best = bound_flags(cand[b, ], bounds),
                  seed_draw = SEED, fold_seeds = paste(SEEDS, collapse = ","), fit_seed = FIT_SEED, stringsAsFactors = FALSE)
  res[[ind]] <- r
  saveRDS(list(candidates = cand, bounds = bounds, tunes = lapply(tn, function(t) t[c("candidates", "cv_by_fold", "fold_assignment", "overfit")]),
               sets = sets, picks = picks, loss_rep = loss_rep, ratio_rep = ratio_rep), file.path(OUT, "detail", paste0(ind, ".rds")))
  utils::write.csv(do.call(rbind, res), file.path(OUT, sprintf("followup_partial_%s.csv", Sys.getenv("TAG", "a"))), row.names = FALSE)
  cat(ind, "done\n"); print(r[, c("indicator", "widening", "pick_row", "new_pick_equals_cv3_pick", "boundary_flags_new_pick")], row.names = FALSE)
}
## the final results csv: Sys.getenv("FINAL") = "1" on the last instance, or concatenate the partial files by hand
if (Sys.getenv("FINAL") == "1") utils::write.csv(do.call(rbind, res), file.path(BASE, "drc_retune_followup.csv"), row.names = FALSE)
