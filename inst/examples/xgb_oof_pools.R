# Out-of-fold residual pools in xgb(): tune once, reuse the tuning folds for the area pool.
# (Shown on the eusilc data of the package; takes a few minutes, so it is not run by the checks.)
library(povmap)
data("eusilcA_pop"); data("eusilcA_smp")
fx <- eqIncome ~ eqsize + cash + self_empl + unempl_ben + age_ben + rent + cap_inv + tax_adj
tune <- xgb_tune(fx, smp_data = eusilcA_smp, domains = "district", transformation = "log",
                 search = "random", n_iter = 20, folds = 10, seed = 1, fit_seed = 123,
                 overfit_ratio = TRUE, keep_oof = TRUE, nthread = 1, verbose = FALSE)
cfg <- xgb_top_configs(tune, rule = "least_overfit", ratio_pool = "domain")   # ranks on the area (domain) pool
est <- xgb(fx, smp_data = eusilcA_smp, pop_data = eusilcA_pop, domains = "district",
           sub_domains = "subDomain", transformation = "log", B = 50, configs = cfg,
           area_pool = "out_of_fold", pool_tune = tune,   # area pool from the tuning folds, no new fits
           sub_pool = "out_of_fold",                      # sub-area pool: fits on sub-area folds
           area_deflate = TRUE, area_sampling_var = NULL) # sampling variance computed from the sample
est$pools$configs[[1]]$area_deflation$factor
# Supplying sampling variances (response scale, squared standard errors):
#   area_sampling_var = c("Wien" = 12.5, "Tirol" = 30.1, ...)   # named by domain label
#   sub_sampling_var  = "my_variance_column"                    # one value per row of smp_data
