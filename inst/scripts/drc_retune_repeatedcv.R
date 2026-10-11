# DRC retune with widened bounds and repeated domain-grouped CV (fold stability), through povmap only.
#
# Candidates: 150 random draws, seed 20261010, bounds as drc_retune.R except lambda 0.01-10000 and min_child_weight 1-1000. Only the
# lambda and min_child_weight values change with the bounds: every parameter is drawn from the same uniform stream, so the other
# parameters' draws are identical (checked in the run log). The draws come out of xgb_tune(search = "random", seed = 20261010), which also
# scores them on the folds of seed 20261010 ("single seed 20261010"). The SAME candidate table is then scored on three further
# domain-grouped 10-fold assignments, seeds 123, 124 and 125 (xgb_tune(grid = candidates, seed = s), fit_seed 123, response scale).
# Analyses, each: paired one-SE set (floor 1, no cap) then least-overfit on the area (domain = territory) pool, ratio_cap NULL:
#   single_20261010, single_123, single_124, single_125   one tune each (its ratios on its own folds)
#   repeated                                              xgb_top_configs(list(t123, t124, t125)): the existing multi-split rule.
#     Loss: mean over the three seeds of each candidate's mean fold loss (= mean over the 30 fold scores).
#     Paired difference from the best: the mean over seeds of the mean fold-wise difference (= mean over the 30 differences).
#     SE: sqrt(mean over seeds of var(fold-wise difference) / 10), the single-split standard error, NOT divided further by the 3 seeds
#     (the three assignments reuse the same sample, so the 30 differences are not independent; this is the package's documented treatment).
#     Ratio: the area-pool out-of-fold / in-sample SD ratio computed on each seed's folds, then averaged over the three seeds.
# Ratios are computed only for the union of the members of the five sets (xgb_overfit_ratio(rows = union)), on each seed's folds; every
# candidate that can be chosen has a ratio for all seeds. Set B is re-scored on every seed's folds as its own one-row tune (not a candidate).
# Environment: INDICATORS, NW (workers), OUT (folder), SEED (20261010), DRC_DATA.
suppressMessages(library(povmap))
DATA <- Sys.getenv("DRC_DATA", "H:/david/SAE/DRC/data")
OUT <- Sys.getenv("OUT", file.path(DATA, "candidate_estimates", "drc_retune_20261010", "repeatedcv"))
N_ITER <- as.integer(Sys.getenv("N_ITER", "150")); NW <- as.integer(Sys.getenv("NW", "5")); SEED <- as.integer(Sys.getenv("SEED", "20261010")); FIT_SEED <- 123L
SEEDS <- c(123L, 124L, 125L)
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
if (nzchar(Sys.getenv("INDICATORS"))) cfg <- cfg[strsplit(Sys.getenv("INDICATORS"), ",")[[1]]]
KEYS <- c("nround", "max_depth", "colsample_bytree", "colsample_bylevel", "colsample_bynode", "subsample", "min_child_weight", "eta", "gamma", "max_delta_step", "lambda", "alpha")
SHOW <- c("nround", "max_depth", "eta", "min_child_weight", "lambda", "subsample", "colsample_bytree")
BOUNDS <- list(eta = c(0.01, 0.3), max_depth = c(2, 10), min_child_weight = c(1, 1000), subsample = c(0.5, 1), colsample_bytree = c(0.3, 1),
               lambda = c(0.01, 10000), nround = c(50, 1500))
BOUNDS_OLD <- BOUNDS; BOUNDS_OLD$min_child_weight <- c(1, 200); BOUNDS_OLD$lambda <- c(0.01, 1000)
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
bound_flags <- function(v) {                     # parameters whose value lies within 5% of the (log or linear) range of an edge of the NEW bounds
  pos <- povmap:::.xgb_bounds_position(BOUNDS, as.list(v))
  f <- pos[pos$at_boundary, ]
  if (!nrow(f)) return("")
  paste(sprintf("%s=%s (%s bound %s)", f$parameter, signif(f$selected, 4), ifelse(f$position < 0.5, "lower", "upper"), ifelse(f$position < 0.5, f$lower, f$upper)), collapse = "; ")
}
res <- list()
for (ind in names(cfg)) {
  tr <- cfg[[ind]][1]; wc <- cfg[[ind]][2]; hc <- cfg[[ind]][3]; hv <- if (is.na(hc)) NULL else hc
  smp <- survey[!is.na(survey[[ind]]), unique(c(ind, X, wc, hv, "ADM2CD_c"))]
  fx <- stats::reformulate(X, response = ind)
  common <- list(fx, smp_data = smp, smp_weights = wc, domains = "ADM2CD_c", transformation = tr, folds = 10, fit_seed = FIT_SEED,
                 nthread = 1, variance_y = hv, score_scale = "response", cpus = NW, verbose = FALSE)
  t0 <- retry(suppressWarnings(do.call(xgb_tune, c(common, list(seed = SEED, search = "random", n_iter = N_ITER, bounds = BOUNDS)))))
  cand <- t0$candidates[, KEYS]; n <- nrow(cand)
  ## the draws other than lambda and min_child_weight are those of the previous bounds
  set.seed(SEED); nc <- length(unique(smp$ADM2CD_c)); invisible(sample(1:10, nc, replace = TRUE))
  d_old <- povmap:::.xgb_sample_bounds(BOUNDS_OLD, N_ITER); same_other <- all(vapply(setdiff(names(d_old), c("lambda", "min_child_weight")),
                                                                             function(p) isTRUE(all.equal(d_old[[p]], cand[[p]])), TRUE))
  tn <- list(); tn[["s20261010"]] <- t0
  for (s in SEEDS) tn[[paste0("s", s)]] <- retry(suppressWarnings(do.call(xgb_tune, c(common, list(seed = s, grid = cand)))))
  rep_tn <- tn[paste0("s", SEEDS)]
  sets <- c(lapply(tn, rows_of, n = n), list(repeated = rows_of(unname(rep_tn), n)))
  U <- sort(unique(unlist(sets)))
  allseeds <- c(SEED, SEEDS)
  for (s in allseeds) {
    nm <- paste0("s", s)
    rt <- retry(xgb_overfit_ratio(fx, smp, cand, smp_weights = wc, domains = "ADM2CD_c", transformation = tr, folds = 10, variance_y = hv,
                                  rows = U, seed = s, fit_seed = FIT_SEED, cpus = NW))
    ov <- data.frame(row = seq_len(n), ratio_total = NA_real_, ratio_domain = NA_real_, ratio_within_domain = NA_real_)
    ov[rt$row, c("ratio_total", "ratio_domain", "ratio_within_domain")] <- rt[, c("ratio_total", "ratio_domain", "ratio_within_domain")]
    tn[[nm]]$overfit <- ov
  }
  rep_tn <- tn[paste0("s", SEEDS)]
  pk <- function(t) attr(xgb_top_configs(t, rule = "least_overfit", ratio_pool = "domain", ratio_cap = NULL), "rows")
  picks <- c(lapply(tn, pk), list(repeated = pk(unname(rep_tn))))
  ## repeated-CV loss and ratios of every candidate
  loss_rep <- rowMeans(vapply(rep_tn, function(t) t$candidates$mse, numeric(n)))
  ratio_rep <- lapply(c("ratio_domain", "ratio_total", "ratio_within_domain"), function(k) rowMeans(vapply(rep_tn, function(t) t$overfit[[k]], numeric(n))))
  names(ratio_rep) <- c("area", "total", "within")
  p <- picks$repeated; b <- which.min(loss_rep)
  ## set B on every seed's folds
  gB <- as.data.frame(setB[[ind]][KEYS])
  tB <- lapply(allseeds, function(s) retry(suppressWarnings(do.call(xgb_tune, c(common, list(seed = s, grid = gB, overfit_ratio = TRUE))))))
  names(tB) <- paste0("s", allseeds)
  tBr <- tB[paste0("s", SEEDS)]
  bm <- function(k) mean(vapply(tBr, function(t) t$overfit[[k]][1], 1))
  r <- data.frame(indicator = ind, transformation = tr, n_candidates = n,
                  set_size_repeated = length(sets$repeated), set_size_s20261010 = length(sets$s20261010), set_size_s123 = length(sets$s123),
                  set_size_s124 = length(sets$s124), set_size_s125 = length(sets$s125),
                  pick_row_repeated = p, cand[p, SHOW],
                  pick_row_s20261010 = picks$s20261010, pick_row_s123 = picks$s123, pick_row_s124 = picks$s124, pick_row_s125 = picks$s125,
                  repeated_equals_s123 = p == picks$s123, repeated_equals_s124 = p == picks$s124, repeated_equals_s125 = p == picks$s125,
                  repeated_equals_s20261010 = p == picks$s20261010,
                  single_picks_all_equal = length(unique(unlist(picks[paste0("s", SEEDS)]))) == 1L,
                  pick_cv_loss_repeated = loss_rep[p], best_row_repeated = b, best_cv_loss_repeated = loss_rep[b],
                  pick_ratio_area = ratio_rep$area[p], pick_ratio_total = ratio_rep$total[p], pick_ratio_within = ratio_rep$within[p],
                  best_ratio_area = ratio_rep$area[b], best_ratio_total = ratio_rep$total[b], best_ratio_within = ratio_rep$within[b],
                  setB_nround = gB$nround, setB_max_depth = gB$max_depth, setB_eta = gB$eta, setB_min_child_weight = gB$min_child_weight,
                  setB_lambda = gB$lambda, setB_subsample = gB$subsample, setB_colsample_bytree = gB$colsample_bytree,
                  setB_cv_loss_repeated = mean(vapply(tBr, function(t) t$candidates$mse[1], 1)),
                  setB_ratio_area = bm("ratio_domain"), setB_ratio_total = bm("ratio_total"), setB_ratio_within = bm("ratio_within_domain"),
                  pick_equals_setB = isTRUE(all(abs(as.numeric(cand[p, SHOW]) - as.numeric(gB[1, SHOW])) <= 1e-9 * pmax(1, abs(as.numeric(gB[1, SHOW]))))),
                  boundary_flags_pick = bound_flags(cand[p, ]), boundary_flags_best = bound_flags(cand[b, ]),
                  other_draws_unchanged_by_bounds = same_other, seed_draw_and_single = SEED, fold_seeds = paste(SEEDS, collapse = ","), fit_seed = FIT_SEED,
                  stringsAsFactors = FALSE)
  res[[ind]] <- r
  saveRDS(list(candidates = cand, tunes = lapply(tn, function(t) t[c("candidates", "cv_by_fold", "fold_assignment", "overfit")]), setB = tB,
               sets = sets, union_rows = U, picks = picks, loss_rep = loss_rep, ratio_rep = ratio_rep), file.path(OUT, "detail", paste0(ind, ".rds")))
  utils::write.csv(do.call(rbind, res), file.path(OUT, sprintf("repeatedcv_partial_%s.csv", Sys.getenv("TAG", "a"))), row.names = FALSE)
  cat(ind, "done\n"); print(r[, c("indicator", "set_size_repeated", "pick_row_repeated", "pick_row_s123", "pick_row_s124", "pick_row_s125", "pick_row_s20261010", "pick_equals_setB", "boundary_flags_pick")], row.names = FALSE)
}
