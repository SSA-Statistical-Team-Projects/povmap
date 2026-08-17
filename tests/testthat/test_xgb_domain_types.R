## Regression test for the domain-type round-trip in xgb().
##
## ?xgb documents `domains` as "can be numeric or a factor". Before the fix,
## point_estim_xgb() rebuilt the domain column from rownames() -- always
## character -- and left_join()ed it straight back onto smp_data, so a numeric
## or integer domain identifier raised:
##   "Can't join `x$<domains>` with `y$<domains>` due to incompatible types."
## The failure was in the point-estimate path, so it hit every call, not only
## bootstrapped ones.

make_toy <- function(domain_type = c("integer", "numeric", "character", "factor"),
                     n_dom = 6L, n_sub = 4L, n_per = 5L, seed = 1L) {
  domain_type <- match.arg(domain_type)
  set.seed(seed)
  pop <- expand.grid(d = seq_len(n_dom), s = seq_len(n_sub))
  pop$sub_id <- paste0("d", pop$d, "s", pop$s)
  pop$x1 <- stats::rnorm(nrow(pop))
  pop$x2 <- stats::runif(nrow(pop))
  pop$popw <- stats::runif(nrow(pop), 5, 50)

  smp <- pop[rep(seq_len(nrow(pop)), each = n_per), ]
  smp$y <- 0.3 * smp$x1 + 0.5 * smp$x2 + stats::rnorm(nrow(smp), sd = 0.1)
  smp$smpw <- stats::runif(nrow(smp), 1, 3)

  cast <- switch(domain_type,
                 integer   = function(v) as.integer(v),
                 numeric   = function(v) as.numeric(v),
                 character = function(v) as.character(v),
                 factor    = function(v) factor(as.character(v)))
  pop$dom <- cast(pop$d); smp$dom <- cast(smp$d)
  list(pop = pop[, c("dom", "sub_id", "x1", "x2", "popw")],
       smp = smp[, c("dom", "sub_id", "x1", "x2", "y", "smpw")])
}

fit_toy <- function(dat) {
  povmap::xgb(fixed = y ~ x1 + x2,
              smp_data = dat$smp, smp_weights = "smpw",
              pop_data = dat$pop, pop_weights = "popw",
              domains = "dom", sub_domains = "sub_id",
              transformation = "no", bootstrap = FALSE,
              nrounds = 5, seed = 1, benchmark_weights = "smpw")
}

test_that("xgb() accepts an integer domain identifier", {
  fit <- fit_toy(make_toy("integer"))
  expect_s3_class(fit, "xgb")
  expect_equal(nrow(fit$ind), 6L)
  expect_false(any(is.na(fit$ind$Mean)))
})

test_that("xgb() accepts a numeric domain identifier", {
  fit <- fit_toy(make_toy("numeric"))
  expect_equal(nrow(fit$ind), 6L)
  expect_false(any(is.na(fit$ind$Mean)))
})

test_that("xgb() accepts a factor domain identifier", {
  fit <- fit_toy(make_toy("factor"))
  expect_equal(nrow(fit$ind), 6L)
  expect_false(any(is.na(fit$ind$Mean)))
})

test_that("domain type does not change the point estimates", {
  # the identifier's storage mode is bookkeeping; the numbers must not move
  a <- fit_toy(make_toy("integer"))$ind
  b <- fit_toy(make_toy("character"))$ind
  a <- a[order(as.character(a$Domain)), ]
  b <- b[order(as.character(b$Domain)), ]
  expect_equal(a$Mean, b$Mean, tolerance = 1e-12)
})
