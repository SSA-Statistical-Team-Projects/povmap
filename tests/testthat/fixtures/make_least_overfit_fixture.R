# Builds least_overfit_drc_fixture.rds from the DRC saved objects (read only). Run by hand.
# For each of the 16 DRC indicators: the candidate rows the DRC selection computed ratios for
# (the members of the one-SE set among the 40 best by fresh score, plus the previous selection and
# the best-scoring configuration), their twelve hyperparameters, fresh-fold score, total ratio and set
# flag, and the hyperparameters of the configuration the DRC selected (xgb_all_tunes_least_overfit).
base <- "H:/david/SAE/DRC/data/candidate_estimates"
KEYS <- c("nround", "max_depth", "colsample_bytree", "colsample_bylevel", "colsample_bynode", "subsample",
          "min_child_weight", "eta", "gamma", "max_delta_step", "lambda", "alpha")
det <- readRDS(file.path(base, "retuned_territory_20261008", "retune_search_detail.rds"))
mem <- read.csv(file.path(base, "least_overfit_20261008", "least_overfit_all_members.csv"))
sel <- read.csv(file.path(base, "least_overfit_20261008", "least_overfit_selection.csv"))
tun <- readRDS(file.path(base, "least_overfit_20261008", "xgb_all_tunes_least_overfit"))
out <- lapply(setNames(names(tun), names(tun)), function(ind) {
  m <- mem[mem$indicator == ind, ]
  cand <- cbind(row = m$row, det[[ind]]$G[m$row, KEYS], fresh_mse = m$fresh_mse,
                ratio_total = m$ratio_total, in_set = m$in_set)
  rownames(cand) <- NULL
  list(cand = cand, expected = as.data.frame(tun[[ind]][KEYS]),
       n_in_set = sel$n_in_set[sel$indicator == ind])
})
saveRDS(out, file.path("tests/testthat/fixtures", "least_overfit_drc_fixture.rds"), version = 2)
