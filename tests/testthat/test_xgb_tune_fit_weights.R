# xgb_tune scores configurations with xgb()'s own case weights (rescaled within domain,
# normalised to mean 1). On raw survey weights, whose leaf sums run into the thousands,
# min_child_weight had no effect on the CV score and lambda almost none.

sim_tune_data <- function() {
  set.seed(42)
  D <- 40; n <- 8
  d <- data.frame(dom = rep(sprintf("d%02d", 1:D), each = n))
  d$x1 <- rnorm(nrow(d)); d$x2 <- rnorm(nrow(d))
  u <- rnorm(D, sd = 0.5)[as.integer(factor(d$dom))]
  d$y <- 1 + d$x1 - 0.5 * d$x2^2 + u + rnorm(nrow(d), sd = 0.5)
  d$w <- runif(nrow(d), 500, 5000)          # survey-scale weights
  d
}

tune_once <- function(d, weight_col = "w") {
  suppressMessages(xgb_tune(y ~ x1 + x2, smp_data = d, smp_weights = weight_col, domains = "dom",
                            folds = 5, nround = 30, max_depth = 3, colsample_bytree = 1,
                            subsample = 1, min_child_weight = c(1, 25), eta = 0.3,
                            lambda = c(0, 25), seed = 7L, cpus = 1, verbose = FALSE))
}

test_that("min_child_weight and lambda produce distinct CV scores", {
  t <- tune_once(sim_tune_data())
  c <- t$candidates
  for (lam in c(0, 25)) {
    s <- c$mse[c$lambda == lam]
    expect_true(abs(diff(s)) > 1e-8 * max(s), info = paste("min_child_weight has no effect at lambda", lam))
  }
  for (mcw in c(1, 25)) {
    s <- c$mse[c$min_child_weight == mcw]
    expect_true(abs(diff(s)) > 1e-8 * max(s), info = paste("lambda has no effect at min_child_weight", mcw))
  }
})

test_that("CV scores do not depend on the scale of the survey weights", {
  d <- sim_tune_data(); d$w1000 <- d$w * 1000
  a <- tune_once(d, "w"); b <- tune_once(d, "w1000")
  expect_equal(a$candidates$mse, b$candidates$mse, tolerance = 1e-8)
  expect_equal(unlist(a[c("min_child_weight", "lambda")]), unlist(b[c("min_child_weight", "lambda")]))
})
