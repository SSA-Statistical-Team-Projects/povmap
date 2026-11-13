# Least-overfit selection, re-run end to end with povmap code only (DRC application).
#
# Reproduces the DRC "set B" selection (DRC_SAE_programs/least_overfit_select.R) from the saved candidate
# tables, using xgb_tune(), xgb_overfit_ratio() and xgb_top_configs(rule = "least_overfit"):
#   1. candidates: the 129 configurations of the territory re-tune (the current tune, the 32-point grid and
#      96 random wider-range draws), saved in retune_search_detail.rds as $G;
#   2. two fresh territory-grouped 10-fold assignments (seeds 901 and 902): xgb_tune(grid = G, seed = s,
#      fit_seed = 123, nthread = 1, score_scale = "response", keep_candidates = TRUE), scoring the weighted
#      territory means on the response scale; model weights are the survey weight rescaled within territory,
#      times variance_y^-0.5 where the indicator has one;
#   3. the paired one-SE set of the two tunes (floor 1, no cap), best 40 by score;
#   4. xgb_overfit_ratio(rows = those 40, seed = 901, fit_seed = 123): ratios in the total, domain and
#      within-domain pools;
#   5. xgb_top_configs(rule = "least_overfit", ratio_pool = "total") on the two tunes: the single
#      configuration with the lowest total ratio.
# Result (DRC, 2026-10-10): 13 of 16 indicators select the identical configuration to set B. access_elec, telpor and
# petitrum differ, because the DRC script took as the set every candidate scoring no worse than the worst member
# of the paired one-SE set (a score prefix; 39, 20 and 28 candidates) while xgb_top_configs searches the paired set
# itself (13, 13 and 17 candidates), and the DRC pick (grid row 104 in each) is outside the paired set. Scores,
# set-prefix sizes and ratios agree to the last digit, so nothing else differs.
# Predictors, province dummies and sample filters are those of the DRC scripts.
#
# Environment variables: DRC_DATA (default H:/david/SAE/DRC/data), OUT (output folder), INDICATORS
# (comma separated, default all 16), NW (workers, default 6), REF (reference tunes to compare with, default
# candidate_estimates/least_overfit_20261008/xgb_all_tunes_least_overfit). Run with povmap and xgboost 3.1.2.1
# on the library path of the workers (R_LIBS).
# Output (OUT): least_overfit_rescore.csv (one row per indicator: selected hyperparameters, the three
# ratios, set size, match with REF) and least_overfit_rescore_detail.rds.
suppressMessages(library(povmap))
DATA <- Sys.getenv("DRC_DATA", "H:/david/SAE/DRC/data")
OUT <- Sys.getenv("OUT", file.path(DATA, "candidate_estimates", "least_overfit_20261008", "e2e_package_rescore"))
NW <- as.integer(Sys.getenv("NW", "6"))
REF <- Sys.getenv("REF", file.path(DATA, "candidate_estimates", "least_overfit_20261008", "xgb_all_tunes_least_overfit"))
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
survey <- readRDS(file.path(DATA, "survey_with_features")); census <- readRDS(file.path(DATA, "census_xgb_agepop"))
X <- colnames(survey)[c(which(colnames(survey) == "population"), which(colnames(survey) == "ntl"):which(colnames(survey) == "FATALITIES"),
                        which(colnames(survey) == "population_territory"), which(colnames(survey) == "ntl_territory"):length(colnames(survey)))]
X <- intersect(setdiff(X, c("popwt", "hhwt", "pop_weight", "wt_5_11", "wt_12_17", "wt_15_64")), names(census))
prov_map <- readRDS(file.path(DATA, "prov_fe_map.rds")); prov_cols <- names(prov_map)
for (pc in prov_cols) survey[[pc]] <- as.integer(as.character(survey$NAM_1) == prov_map[[pc]])
X <- c(X, prov_cols)
det <- readRDS(file.path(DATA, "candidate_estimates", "retuned_territory_20261008", "retune_search_detail.rds"))
cfg <- list(nivie = c("arcsin", "popwt", NA), theil = c("log", "popwt", NA), pcexp = c("log", "popwt", NA), access_elec = c("arcsin", "popwt", NA), access_water = c("arcsin", "popwt", NA),
            lfp = c("arcsin", "wt_15_64", "lfp_var"), school_5_11 = c("arcsin", "wt_5_11", "school_5_11_var"), school_12_17 = c("arcsin", "wt_12_17", "school_12_17_var"),
            telpor = c("arcsin", "popwt", NA), grosrum = c("poisson", "hhwt", NA), petitrum = c("poisson", "hhwt", NA), porc = c("poisson", "hhwt", NA), lapin = c("poisson", "hhwt", NA),
            volail = c("poisson", "hhwt", NA), car = c("arcsin", "popwt", NA), tv = c("arcsin", "popwt", NA))
if (nzchar(Sys.getenv("INDICATORS"))) cfg <- cfg[strsplit(Sys.getenv("INDICATORS"), ",")[[1]]]
KEYS <- c("nround", "max_depth", "colsample_bytree", "colsample_bylevel", "colsample_bynode", "subsample", "min_child_weight", "eta", "gamma", "max_delta_step", "lambda", "alpha")
ref <- if (file.exists(REF)) readRDS(REF) else NULL
NCAP <- 40L
## starting a PSOCK cluster occasionally fails on busy Windows hosts ("port cannot be opened"): retry
retry <- function(expr, tries = 6) {
  e <- substitute(expr); pf <- parent.frame()
  for (i in seq_len(tries)) {
    r <- tryCatch(eval(e, pf), error = function(err) if (grepl("serverSocket|port", conditionMessage(err)) && i < tries) NULL else stop(err))
    if (!is.null(r)) return(r)
    Sys.sleep(5)
  }
}
res <- list(); detail <- list()
for (ind in names(cfg)) {
  tr <- cfg[[ind]][1]; wc <- cfg[[ind]][2]; hc <- cfg[[ind]][3]; hv <- if (is.na(hc)) NULL else hc
  G <- det[[ind]]$G[, KEYS]
  smp <- survey[!is.na(survey[[ind]]), unique(c(ind, X, wc, hv, "ADM2CD_c"))]
  fx <- stats::reformulate(X, response = ind)
  tn <- lapply(c(901, 902), function(s) retry(suppressWarnings(xgb_tune(
    fx, smp_data = smp, smp_weights = wc, domains = "ADM2CD_c", transformation = tr, folds = 10, grid = G,
    seed = s, fit_seed = 123, nthread = 1, variance_y = hv, score_scale = "response", cpus = NW, verbose = FALSE))))
  full <- xgb_top_configs(tn, rule = "one_se_paired", floor = 1, cap = nrow(G))
  rows_set <- {
    sc <- (rowMeans(tn[[1]]$cv_by_fold) + rowMeans(tn[[2]]$cv_by_fold)) / 2
    r <- which(sc <= max(full$cv_score) + 1e-15); r[order(sc[r])]
  }
  ## the paired one-SE set itself (all members of the rule, ratio-free): a diagnostic. The DRC script took as the set
  ## every candidate scoring no worse than the worst member (a score prefix), which can contain candidates the
  ## paired rule excludes.
  paired_rows <- attr(xgb_top_configs(tn, rule = "one_se_ratio", ratio = rep(1, nrow(G)), max_ratio = 2, cap = nrow(G),
                                      ratio_cap = NULL), "rows")
  used <- rows_set[seq_len(min(NCAP, length(rows_set)))]
  rt <- retry(xgb_overfit_ratio(fx, smp, G, smp_weights = wc, domains = "ADM2CD_c", transformation = tr, folds = 10,
                          variance_y = hv, rows = used, seed = 901, fit_seed = 123, cpus = NW))
  pick <- xgb_top_configs(tn, rule = "least_overfit", ratio = rt, ratio_pool = "total", ratio_cap = NCAP)
  prow <- attr(pick, "rows")
  got <- G[prow, KEYS]
  r <- data.frame(indicator = ind, n_in_set = length(rows_set), n_paired_set = length(paired_rows), n_ratioed = length(used), grid_row = prow, got,
                  ratio_total = rt$ratio_total[rt$row == prow], ratio_domain = rt$ratio_domain[rt$row == prow],
                  ratio_within_domain = rt$ratio_within_domain[rt$row == prow])
  if (!is.null(ref)) {
    r$matches_ref <- isTRUE(all(as.numeric(got[1, ]) == as.numeric(unlist(ref[[ind]][KEYS]))))
    r$max_abs_diff_ref <- max(abs(as.numeric(got[1, ]) - as.numeric(unlist(ref[[ind]][KEYS]))))
  }
  res[[ind]] <- r; detail[[ind]] <- list(set_rows = rows_set, paired_rows = paired_rows, ratios = rt, pick = pick,
                                          fresh_mse = (rowMeans(tn[[1]]$cv_by_fold) + rowMeans(tn[[2]]$cv_by_fold)) / 2)
  cat(ind, "done\n"); print(r, row.names = FALSE)
  utils::write.csv(do.call(rbind, res), file.path(OUT, "least_overfit_rescore.csv"), row.names = FALSE)
  saveRDS(detail, file.path(OUT, "least_overfit_rescore_detail.rds"))
}
