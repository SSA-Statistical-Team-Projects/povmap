# DRC retune of the 16 indicators under the shared XGBoost convention (Stage 1), through povmap::xgb_tune only.
#
# For each indicator:
#   xgb_tune(search = "random", n_iter = 150, bounds = convention bounds, folds = 10 grouped by ADM2CD_c, seed = SEED,
#            fit_seed = FIT_SEED, nthread = 1, score_scale = "response", variance_y as in the DRC models,
#            overfit_ratio = TRUE)                     # fold assignment is saved in the result ($fold_assignment)
#   xgb_top_configs(rule = "one_se_paired", floor = 1, cap = all)      -> candidate set (size recorded)
#   xgb_top_configs(rule = "least_overfit", ratio_pool = "domain", ratio_cap = NULL) -> the pick (area = territory pool)
# The set B configuration (DRC, 2026-10-08) is re-scored on the SAME folds (xgb_tune(grid = set B row), same seed, same
# fit seed) so its CV loss and overfit ratios are comparable with the new pick's.
# Environment: DRC_DATA, OUT (default the shared drc_retune_20261010 folder), INDICATORS, NW (workers, default 6),
# SEED (default 20261010), FIT_SEED (123). Needs povmap (feature/least-overfit-selection) and xgboost 3.1.2.1 on R_LIBS.
suppressMessages(library(povmap))
DATA <- Sys.getenv("DRC_DATA", "H:/david/SAE/DRC/data")
OUT <- Sys.getenv("OUT", file.path(DATA, "candidate_estimates", "drc_retune_20261010"))
N_ITER <- as.integer(Sys.getenv("N_ITER", "150")); NW <- as.integer(Sys.getenv("NW", "6")); SEED <- as.integer(Sys.getenv("SEED", "20261010")); FIT_SEED <- as.integer(Sys.getenv("FIT_SEED", "123"))
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
## convention search space: package default bounds (Colombia's), with lambda 0.01-1000 and min_child_weight 1-200
BOUNDS <- list(eta = c(0.01, 0.3), max_depth = c(2, 10), min_child_weight = c(1, 200), subsample = c(0.5, 1), colsample_bytree = c(0.3, 1),
               lambda = c(0.01, 1000), nround = c(50, 1500))
setB <- readRDS(file.path(DATA, "candidate_estimates", "least_overfit_20261008", "xgb_all_tunes_least_overfit"))
retry <- function(expr, tries = 6) {            # starting a PSOCK cluster occasionally fails on busy Windows hosts
  e <- substitute(expr); pf <- parent.frame()
  for (i in seq_len(tries)) {
    r <- tryCatch(eval(e, pf), error = function(err) if (grepl("serverSocket|port", conditionMessage(err)) && i < tries) NULL else stop(err))
    if (!is.null(r)) return(r)
    Sys.sleep(5)
  }
}
## set B outside the search space? (stop if so)
for (ind in names(cfg)) {
  b <- unlist(setB[[ind]][c("nround", "max_depth", "eta", "min_child_weight", "lambda", "subsample", "colsample_bytree")])
  lo <- vapply(BOUNDS, `[`, 1, FUN.VALUE = 1)[names(b)]; hi <- vapply(BOUNDS, `[`, 2, FUN.VALUE = 1)[names(b)]
  if (any(b < lo | b > hi)) stop("set B for ", ind, " lies outside the search space: ", paste(names(b)[b < lo | b > hi], collapse = ", "))
}
res <- list()
for (ind in names(cfg)) {
  tr <- cfg[[ind]][1]; wc <- cfg[[ind]][2]; hc <- cfg[[ind]][3]; hv <- if (is.na(hc)) NULL else hc
  smp <- survey[!is.na(survey[[ind]]), unique(c(ind, X, wc, hv, "ADM2CD_c"))]
  fx <- stats::reformulate(X, response = ind)
  common <- list(fx, smp_data = smp, smp_weights = wc, domains = "ADM2CD_c", transformation = tr, folds = 10, seed = SEED, fit_seed = FIT_SEED,
                 nthread = 1, variance_y = hv, score_scale = "response", overfit_ratio = TRUE, cpus = NW, verbose = FALSE)
  edge_warn <- character()
  tn <- retry(withCallingHandlers(do.call(xgb_tune, c(common, list(search = "random", n_iter = N_ITER, bounds = BOUNDS))),
                                  warning = function(w) { edge_warn <<- c(edge_warn, conditionMessage(w)); invokeRestart("muffleWarning") }))
  n <- nrow(tn$candidates)
  full <- xgb_top_configs(tn, rule = "one_se_paired", floor = 1, cap = n)
  pick <- xgb_top_configs(tn, rule = "least_overfit", ratio_pool = "domain", ratio_cap = NULL)
  prow <- attr(pick, "rows"); brow <- which.min(tn$candidates$mse)
  ## set B on the same folds, same fit seed
  gB <- as.data.frame(setB[[ind]][KEYS])
  tB <- retry(suppressWarnings(do.call(xgb_tune, c(common, list(grid = gB)))))
  ov <- tn$overfit
  r <- data.frame(indicator = ind, transformation = tr, n_candidates = n, set_size = nrow(full),
                  pick_row = prow, pick[1, c("nrounds", "max_depth", "eta", "min_child_weight", "lambda", "subsample", "colsample_bytree",
                                              "colsample_bylevel", "colsample_bynode", "gamma", "max_delta_step", "alpha")],
                  pick_cv_loss = tn$candidates$mse[prow], best_row = brow, best_cv_loss = tn$candidates$mse[brow],
                  pick_ratio_total = ov$ratio_total[prow], pick_ratio_area = ov$ratio_domain[prow], pick_ratio_within = ov$ratio_within_domain[prow],
                  best_ratio_total = ov$ratio_total[brow], best_ratio_area = ov$ratio_domain[brow], best_ratio_within = ov$ratio_within_domain[brow],
                  min_area_ratio_in_set = min(ov$ratio_domain[attr(xgb_top_configs(tn, rule = "one_se_ratio", ratio = rep(1, n), max_ratio = 2, cap = n, ratio_cap = NULL), "rows")]),
                  setB_nrounds = gB$nround, setB_max_depth = gB$max_depth, setB_eta = gB$eta, setB_min_child_weight = gB$min_child_weight,
                  setB_lambda = gB$lambda, setB_subsample = gB$subsample, setB_colsample_bytree = gB$colsample_bytree,
                  setB_cv_loss_same_folds = tB$candidates$mse[1], setB_ratio_total = tB$overfit$ratio_total[1],
                  setB_ratio_area = tB$overfit$ratio_domain[1], setB_ratio_within = tB$overfit$ratio_within_domain[1],
                  same_as_setB = isTRUE(all(as.numeric(pick[1, c("nrounds", "max_depth", "eta", "min_child_weight", "lambda", "subsample", "colsample_bytree")]) ==
                                            as.numeric(gB[1, c("nround", "max_depth", "eta", "min_child_weight", "lambda", "subsample", "colsample_bytree")]))),
                  boundary_warning = if (length(edge_warn)) edge_warn[1] else "", seed = SEED, fit_seed = FIT_SEED)
  res[[ind]] <- r
  saveRDS(list(tune = tn, setB_tune = tB, set_rows = attr(xgb_top_configs(tn, rule = "one_se_ratio", ratio = rep(1, n), max_ratio = 2, cap = n, ratio_cap = NULL), "rows"),
               pick = pick, fold_assignment = tn$fold_assignment), file.path(OUT, "detail", paste0(ind, ".rds")))
  utils::write.csv(do.call(rbind, res), file.path(OUT, "drc_retune_selections_partial.csv"), row.names = FALSE)
  cat(ind, "done\n"); print(r[, c("indicator", "set_size", "pick_row", "same_as_setB", "pick_ratio_area", "setB_ratio_area")], row.names = FALSE)
}
utils::write.csv(do.call(rbind, res), file.path(OUT, "drc_retune_selections.csv"), row.names = FALSE)
