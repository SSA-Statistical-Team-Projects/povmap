# Decomposition of the move from DRC set B to the convention pick (Stage 1), by a ladder of single changes on the saved
# 129-candidate pool of the DRC retune (retune_search_detail.rds $G), all through povmap functions.
#   step 0  set B               : DRC scripts (prefix set, total ratio, folds 901/902)            -> xgb_all_tunes_least_overfit
#   step 1  candidate set       : R1 total   = two fresh fold seeds (901, 902), TRUE paired one-SE set, rank on the total pool
#   step 2  ranking pool        : R1 area    = same, rank on the area (domain = territory) pool
#   step 3  fold protocol       : R2 area    = convention protocol on the same pool: one run, seed SEED, ratios on the search folds
#   step 4  search space        : the convention pick from the 150-draw wide-bounds search (drc_retune.R)
# Scoring scale is response in every step (set B was already scored on the response scale), so that step has no move.
# The poisson fix is isolated for the five count indicators by drc_retune.R with TR_OVERRIDE=no (everything else as the convention run).
# The ladder is sequential: a move is credited to the first step at which the pick changes, and the order matters.
# Output (OUT): drc_ladder_picks.csv (the pick at each step, grid rows and hyperparameters) and ladder_detail.rds.
suppressMessages(library(povmap))
DATA <- Sys.getenv("DRC_DATA", "H:/david/SAE/DRC/data")
OUT <- Sys.getenv("OUT", file.path(DATA, "candidate_estimates", "drc_retune_20261010"))
NW <- as.integer(Sys.getenv("NW", "6")); SEED <- as.integer(Sys.getenv("SEED", "20261010")); FIT_SEED <- 123L
survey <- readRDS(file.path(DATA, "survey_with_features")); census <- readRDS(file.path(DATA, "census_xgb_agepop"))
X <- colnames(survey)[c(which(colnames(survey) == "population"), which(colnames(survey) == "ntl"):which(colnames(survey) == "FATALITIES"),
                        which(colnames(survey) == "population_territory"), which(colnames(survey) == "ntl_territory"):length(colnames(survey)))]
X <- intersect(setdiff(X, c("popwt", "hhwt", "pop_weight", "wt_5_11", "wt_12_17", "wt_15_64")), names(census))
prov_map <- readRDS(file.path(DATA, "prov_fe_map.rds")); prov_cols <- names(prov_map)
for (pc in prov_cols) survey[[pc]] <- as.integer(as.character(survey$NAM_1) == prov_map[[pc]])
X <- c(X, prov_cols)
det <- readRDS(file.path(DATA, "candidate_estimates", "retuned_territory_20261008", "retune_search_detail.rds"))
setB <- readRDS(file.path(DATA, "candidate_estimates", "least_overfit_20261008", "xgb_all_tunes_least_overfit"))
cfg <- list(nivie = c("arcsin", "popwt", NA), theil = c("log", "popwt", NA), pcexp = c("log", "popwt", NA), access_elec = c("arcsin", "popwt", NA), access_water = c("arcsin", "popwt", NA),
            lfp = c("arcsin", "wt_15_64", "lfp_var"), school_5_11 = c("arcsin", "wt_5_11", "school_5_11_var"), school_12_17 = c("arcsin", "wt_12_17", "school_12_17_var"),
            telpor = c("arcsin", "popwt", NA), grosrum = c("poisson", "hhwt", NA), petitrum = c("poisson", "hhwt", NA), porc = c("poisson", "hhwt", NA), lapin = c("poisson", "hhwt", NA),
            volail = c("poisson", "hhwt", NA), car = c("arcsin", "popwt", NA), tv = c("arcsin", "popwt", NA))
if (nzchar(Sys.getenv("INDICATORS"))) cfg <- cfg[strsplit(Sys.getenv("INDICATORS"), ",")[[1]]]
KEYS <- c("nround", "max_depth", "colsample_bytree", "colsample_bylevel", "colsample_bynode", "subsample", "min_child_weight", "eta", "gamma", "max_delta_step", "lambda", "alpha")
retry <- function(expr, tries = 6) {
  e <- substitute(expr); pf <- parent.frame()
  for (i in seq_len(tries)) {
    r <- tryCatch(eval(e, pf), error = function(err) if (grepl("serverSocket|port", conditionMessage(err)) && i < tries) NULL else stop(err))
    if (!is.null(r)) return(r)
    Sys.sleep(5)
  }
}
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
rows_of <- function(tn, n) attr(xgb_top_configs(tn, rule = "one_se_ratio", ratio = rep(1, n), max_ratio = 2, cap = n, ratio_cap = NULL), "rows")
out <- list(); detail <- list()
for (ind in names(cfg)) {
  tr <- cfg[[ind]][1]; wc <- cfg[[ind]][2]; hc <- cfg[[ind]][3]; hv <- if (is.na(hc)) NULL else hc
  G <- det[[ind]]$G[, KEYS]; n <- nrow(G)
  smp <- survey[!is.na(survey[[ind]]), unique(c(ind, X, wc, hv, "ADM2CD_c"))]
  fx <- stats::reformulate(X, response = ind)
  base <- list(fx, smp_data = smp, smp_weights = wc, domains = "ADM2CD_c", transformation = tr, folds = 10, fit_seed = FIT_SEED, nthread = 1,
               variance_y = hv, score_scale = "response", grid = G, cpus = NW, verbose = FALSE)
  ## R1: two fresh fold assignments, true paired set, ratios on seed 901
  tn <- lapply(c(901, 902), function(s) retry(suppressWarnings(do.call(xgb_tune, c(base, list(seed = s))))))
  prs <- rows_of(tn, n)
  rt <- retry(xgb_overfit_ratio(fx, smp, G, smp_weights = wc, domains = "ADM2CD_c", transformation = tr, folds = 10, variance_y = hv,
                                rows = prs, seed = 901, fit_seed = FIT_SEED, cpus = NW))
  p1t <- attr(xgb_top_configs(tn, rule = "least_overfit", ratio = rt, ratio_pool = "total", ratio_cap = NULL), "rows")
  p1a <- attr(xgb_top_configs(tn, rule = "least_overfit", ratio = rt, ratio_pool = "domain", ratio_cap = NULL), "rows")
  ## R2: convention protocol on the same pool
  t2 <- retry(suppressWarnings(do.call(xgb_tune, c(base, list(seed = SEED, overfit_ratio = TRUE)))))
  p2a <- attr(xgb_top_configs(t2, rule = "least_overfit", ratio_pool = "domain", ratio_cap = NULL), "rows")
  ## set B row in the pool (matching hyperparameters)
  b <- as.numeric(unlist(setB[[ind]][KEYS])); rB <- which(apply(G, 1, function(r) all(as.numeric(r) == b)))[1]
  one <- function(step, row) data.frame(indicator = ind, step = step, grid_row = row, G[row, c("nround", "max_depth", "eta", "min_child_weight", "lambda", "subsample", "colsample_bytree")])
  out[[ind]] <- rbind(one("0_setB", rB), one("1_paired_set_total_pool", p1t), one("2_area_pool", p1a), one("3_convention_folds_area_pool", p2a))
  detail[[ind]] <- list(paired_rows_R1 = prs, ratios_R1 = rt, pick = c(rB = rB, p1t = p1t, p1a = p1a, p2a = p2a), setB_row_in_pool = rB,
                        set_rows_R2 = rows_of(t2, n), overfit_R2 = t2$overfit, mse_R2 = t2$candidates$mse)
  cat(ind, "done: setB row", rB, "R1 total", p1t, "R1 area", p1a, "R2 area", p2a, "paired set", length(prs), "\n")
  utils::write.csv(do.call(rbind, out), file.path(OUT, "drc_ladder_picks.csv"), row.names = FALSE)
  saveRDS(detail, file.path(OUT, "ladder_detail.rds"))
}
