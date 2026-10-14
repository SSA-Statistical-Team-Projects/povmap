test_that("ebp with L = 50 runs on eusilcA and returns every indicator column", {
  data("eusilcA_pop", package = "povmap"); data("eusilcA_smp", package = "povmap")
  fx <- eqIncome ~ gender + eqsize + cash + self_empl + unempl_ben + age_ben + surv_ben +
    sick_ben + dis_ben + rent + fam_allow + house_allow + cap_inv + tax_adj
  m <- suppressMessages(ebp(fixed = fx, pop_data = eusilcA_pop, pop_domains = "district",
                            smp_data = eusilcA_smp, smp_domains = "district",
                            L = 50, MSE = FALSE, na.rm = TRUE, cpus = 1))
  expect_s3_class(m, "ebp")
  expect_equal(ncol(m$ind) - 1, length(m$framework$indicator_names))
  expect_true(all(c("Mean", "Minimum", "Maximum", "Quantile_5", "Median", "Quantile_95") %in% names(m$ind)))
  expect_false(anyNA(m$ind$Mean))
})

test_that("a subset of indicators including Quantiles also runs", {
  data("eusilcA_pop", package = "povmap"); data("eusilcA_smp", package = "povmap")
  fx <- eqIncome ~ gender + eqsize + cash + self_empl + unempl_ben + age_ben + surv_ben +
    sick_ben + dis_ben + rent + fam_allow + house_allow + cap_inv + tax_adj
  m <- suppressMessages(ebp(fixed = fx, pop_data = eusilcA_pop, pop_domains = "district",
                            smp_data = eusilcA_smp, smp_domains = "district",
                            L = 20, MSE = FALSE, na.rm = TRUE, cpus = 1,
                            indicators = c("Mean", "Quantiles")))
  expect_equal(names(m$ind), c("Domain", "Mean", paste0("Quantile_", c(5, 10, 25)), "Median",
                               paste0("Quantile_", c(75, 90, 95))))
})
