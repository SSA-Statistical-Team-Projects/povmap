# Reproduce the Stage 1b scratch-hook result (holdout_coverage_20261010, oof_lib.R + job_oof.R) with the package's own
# out-of-fold pools: DRC access_water, CV seed 123, fold 1, B = 100. The scratch job built the pools outside xgb() and injected
# them with a text-substituted copy of xgb(); here xgb(area_pool = "out_of_fold", sub_pool = "out_of_fold_domain") does the same.
# Compares Mean, Mean_bench, Lower, Upper, Lower_bench, Upper_bench of the held-out territories with the saved job result.
# Environment: IND (default access_water), FOLD (1), OUT (comparison csv). Needs povmap + xgboost 3.1.2.1 on R_LIBS.
suppressMessages({library(povmap); library(dplyr); library(haven)})
stopifnot(as.character(packageVersion("xgboost")) == "3.1.2.1")
IND <- Sys.getenv("IND", "access_water"); FOLD <- as.integer(Sys.getenv("FOLD", "1")); CVSEED <- 123L; B <- 100L
DRC_BASE <- "H:/david/SAE/DRC"; DATA_DIR <- file.path(DRC_BASE, "data")
HOLD <- file.path(DATA_DIR, "candidate_estimates", "holdout_coverage_20261010")
OUT <- Sys.getenv("OUT", file.path(DATA_DIR, "candidate_estimates", "drc_retune_20261010", sprintf("oof_pool_repro_%s_f%02d.csv", IND, FOLD)))
survey <- readRDS(file.path(DATA_DIR, "survey_with_features")); census_xgb <- readRDS(file.path(DATA_DIR, "census_xgb_agepop"))
indicators <- c("nivie","theil","pcexp","access_elec","access_water","lfp","school_5_11","school_12_17","telpor","grosrum","petitrum","porc","lapin","volail","car","tv")
livestock <- c("grosrum","petitrum","porc","lapin","volail")
smp_weight <- setNames(rep("popwt", 16), indicators); smp_weight[livestock] <- "hhwt"
smp_weight[c("school_5_11","school_12_17","lfp")] <- c("wt_5_11","wt_12_17","wt_15_64")
pop_weight <- setNames(rep("pop_weight", 16), indicators)
pop_weight[c("school_5_11","school_12_17","lfp")] <- c("population_5_11","population_12_17","population_15_64")
transformations <- setNames(c("arcsin","log","log","arcsin","arcsin","arcsin","arcsin","arcsin","arcsin","poisson","poisson","poisson","poisson","poisson","arcsin","arcsin"), indicators)
benchmark_types <- setNames(c("logit_raking","ratio","ratio","logit_raking","logit_raking","logit_raking","logit_raking","logit_raking","logit_raking","ratio","ratio","ratio","ratio","ratio","logit_raking","logit_raking"), indicators)
Xvars_xgb <- colnames(survey)[c(which(colnames(survey)=="population"), which(colnames(survey)=="ntl"):which(colnames(survey)=="FATALITIES"),
                                which(colnames(survey)=="population_territory"), which(colnames(survey)=="ntl_territory"):length(colnames(survey)))]
Xvars_xgb <- setdiff(Xvars_xgb, c("popwt","hhwt","pop_weight", unname(smp_weight))); Xvars_xgb <- intersect(Xvars_xgb, names(census_xgb))
prov_map <- readRDS(file.path(DATA_DIR, "prov_fe_map.rds")); prov_cols <- names(prov_map)
add_prov_fe <- function(df) { for (pc in prov_cols) df[[pc]] <- as.integer(as.character(df$NAM_1) == prov_map[[pc]]); df }
survey <- add_prov_fe(survey); census_xgb <- add_prov_fe(census_xgb)
het_var <- c(lfp="lfp_var", school_5_11="school_5_11_var", school_12_17="school_12_17_var")
prov_direct <- haven::read_dta(file.path(DRC_BASE, "DRC_SAE_programs/reproducibility package/storage/direct_estimates_province.dta"))
bench_target <- function(ind) { d <- prov_direct[prov_direct$indicator == ind, ]; setNames(d$direct, as.character(d$NAM_1)) }
bench_target_se <- function(ind) { d <- prov_direct[prov_direct$indicator == ind, ]; setNames(d$direct_SD_bench, as.character(d$NAM_1)) }
all_tunes <- readRDS(file.path(DATA_DIR, "candidate_estimates/least_overfit_20261008/xgb_all_tunes_least_overfit"))
Xv <- c(Xvars_xgb, prov_cols); variance_y <- if (IND %in% names(het_var)) het_var[[IND]] else NULL
smp <- survey[!is.na(survey[[IND]]), unique(c(IND, Xv, smp_weight[[IND]], variance_y, "ADM2CD_c", "poly_id", "NAM_1"))]
set.seed(CVSEED); ud <- unique(smp[, "ADM2CD_c"]); fa <- data.frame(Domain = ud, fold = sample(rep(1:10, length.out = length(ud))), stringsAsFactors = FALSE)
hold <- fa$Domain[fa$fold == FOLD]
train <- smp[!(smp[, "ADM2CD_c"] %in% hold), ]
tune <- all_tunes[[IND]]
a <- list(fixed = as.formula(paste0(IND, paste("~", paste(Xv, collapse = " + ")))), smp_data = train, smp_weights = smp_weight[[IND]], pop_data = census_xgb,
  pop_weights = pop_weight[[IND]], domains = "ADM2CD_c", sub_domains = "poly_id", transformation = transformations[[IND]],
  nrounds = tune$nround, max_depth = tune$max_depth, colsample_bytree = tune$colsample_bytree, colsample_bylevel = tune$colsample_bylevel,
  colsample_bynode = tune$colsample_bynode, subsample = tune$subsample, min_child_weight = tune$min_child_weight, eta = tune$eta,
  gamma = tune$gamma, max_delta_step = tune$max_delta_step, lambda = tune$lambda, alpha = tune$alpha, na.rm = TRUE, weightedBS = TRUE,
  benchmark = bench_target(IND), benchmark_level = "NAM_1", benchmark_type = benchmark_types[[IND]],
  variance_y = variance_y, perturb_benchmark = TRUE, benchmark_target_se = bench_target_se(IND),
  B = B, seed = 123, cpus = 1, interval_method = "basic", interval_scale = "natural", bootstrap = TRUE)
est <- do.call(xgb, c(a, list(area_pool = "out_of_fold", sub_pool = "out_of_fold_domain", pool_folds = 10, pool_seed = 123)))
if (transformations[[IND]] == "log") for (cn in grep("^(Lower|Upper)", names(est$CI), value = TRUE)) est$CI[[cn]] <- pmax(est$CI[[cn]], 0)
h <- as.character(est$ind$Domain) %in% as.character(hold)
res <- data.frame(Domain = as.character(est$ind$Domain)[h], Mean = est$ind$Mean[h], Mean_bench = est$ind$Mean_bench[h], Lower = est$CI$Lower[h],
                  Upper = est$CI$Upper[h], Lower_bench = est$CI$Lower_bench[h], Upper_bench = est$CI$Upper_bench[h])
ref <- readRDS(file.path(HOLD, "jobs", sprintf("oof_%s_cv%d_f%02d_B%d.rds", IND, CVSEED, FOLD, B)))$res
ins <- readRDS(file.path(HOLD, "jobs", sprintf("ins_%s_cv%d_f%02d_B%d.rds", IND, CVSEED, FOLD, B)))$res
m <- merge(res, ref, by = "Domain", suffixes = c(".pkg", ".scratch")); m <- merge(m, ins[, c("Domain", "Mean", "Mean_bench")], by = "Domain")
stopifnot(nrow(m) == nrow(ref))
cols <- c("Mean", "Mean_bench", "Lower", "Upper", "Lower_bench", "Upper_bench")
cmp <- data.frame(quantity = cols, max_abs_diff = vapply(cols, function(k) max(abs(m[[paste0(k, ".pkg")]] - m[[paste0(k, ".scratch")]])), 1),
                  n_territories = nrow(m))
cmp$max_abs_diff_vs_insample_point <- c(max(abs(m$Mean.pkg - m$Mean)), max(abs(m$Mean_bench.pkg - m$Mean_bench)), NA, NA, NA, NA)
cmp$max_halfwidth_diff <- c(NA, NA, max(abs((m$Upper.pkg - m$Lower.pkg) - (m$Upper.scratch - m$Lower.scratch)) / 2), NA, NA, NA)
dir.create(dirname(OUT), showWarnings = FALSE, recursive = TRUE); utils::write.csv(cmp, OUT, row.names = FALSE); print(cmp)
