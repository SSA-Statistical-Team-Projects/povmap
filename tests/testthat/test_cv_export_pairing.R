# CV / precision columns are paired with point estimates by name, not position.
# Since c552f02, xgb's $ind carries Mean_agg, which has no variance in $var; the
# benchmarked variance is stored as $var$Var_bench and belongs to $ind$Mean_bench.

xgb_like <- function(with_agg = TRUE) {
  ind <- data.frame(Domain = c("a", "b", "c"), Mean = c(0.2, 0.4, 0.5))
  if (with_agg) ind$Mean_agg <- c(0.21, 0.41, 0.52)
  ind$Mean_bench <- c(0.25, 0.35, 0.55)
  # var rows deliberately in a different order: rows are matched on Domain
  var <- data.frame(Domain = c("c", "a", "b"), Mean = c(0.0025, 0.0004, 0.0016),
                    Var_bench = c(0.0036, 0.0009, 0.0025))
  structure(list(ind = ind, var = var), class = c("xgb", "povmap"))
}

test_that("xgb CV pairs Mean with var$Mean and Mean_bench with var$Var_bench", {
  obj <- xgb_like(with_agg = TRUE)
  out <- estimators(obj, indicator = "all", var = TRUE, CV = TRUE)$ind
  expect_identical(names(out), c("Domain", "Mean", "Mean_agg", "Mean_bench",
                                 "Mean_Var", "Mean_bench_Var", "Mean_CV", "Mean_bench_CV"))
  v <- obj$var[match(out$Domain, obj$var$Domain), ]
  expect_equal(out$Mean_Var, v$Mean)
  expect_equal(out$Mean_bench_Var, v$Var_bench)
  expect_equal(out$Mean_CV, sqrt(v$Mean) / obj$ind$Mean)
  expect_equal(out$Mean_bench_CV, sqrt(v$Var_bench) / obj$ind$Mean_bench)
  expect_equal(out$Mean_agg, obj$ind$Mean_agg)
})

test_that("xgb objects without Mean_agg export as before", {
  obj <- xgb_like(with_agg = FALSE)
  out <- estimators(obj, indicator = "all", var = TRUE, CV = TRUE)$ind
  expect_identical(names(out), c("Domain", "Mean", "Mean_bench", "Mean_Var",
                                 "Mean_bench_Var", "Mean_CV", "Mean_bench_CV"))
  v <- obj$var[match(out$Domain, obj$var$Domain), ]
  expect_equal(out$Mean_bench_CV, sqrt(v$Var_bench) / obj$ind$Mean_bench)
})

test_that("selecting Mean_bench returns its own variance and CV", {
  out <- estimators(xgb_like(), indicator = "Mean_bench", var = TRUE, CV = TRUE)$ind
  expect_identical(names(out), c("Domain", "Mean_bench", "Mean_bench_Var", "Mean_bench_CV"))
  expect_equal(out$Mean_bench_Var, c(0.0009, 0.0025, 0.0036))
})

test_that("mse_emdi (split export path) uses the same pairing", {
  obj <- xgb_like()
  obj$MSE <- obj$var
  p <- povmap:::mse_emdi(obj, indicator = "all", CV = TRUE)
  expect_identical(names(p$ind), c("Domain", "Mean", "Mean_bench"))
  expect_equal(p$ind_cv$Mean_bench, sqrt(obj$var$Var_bench) /
                 obj$ind$Mean_bench[match(obj$var$Domain, obj$ind$Domain)])
})

test_that("models whose MSE columns share the point names are unchanged", {
  obj <- structure(list(
    ind = data.frame(Domain = c("a", "b"), Mean = c(10, 20), Head_Count = c(0.1, 0.2)),
    MSE = data.frame(Domain = c("a", "b"), Mean = c(4, 9), Head_Count = c(0.0001, 0.0004))),
    class = c("ebp", "povmap"))
  out <- estimators(obj, indicator = "all", MSE = TRUE, CV = TRUE)$ind
  expect_identical(names(out), c("Domain", "Mean", "Head_Count", "Mean_MSE",
                                 "Head_Count_MSE", "Mean_CV", "Head_Count_CV"))
  expect_equal(out$Mean_CV, c(0.2, 0.15))
  expect_equal(out$Head_Count_CV, c(0.1, 0.1))
})

test_that("a precision column with no point estimate is an error", {
  obj <- xgb_like()
  obj$var$Extra <- 1
  expect_error(estimators(obj, indicator = "all", var = TRUE, CV = TRUE), "Extra")
})
