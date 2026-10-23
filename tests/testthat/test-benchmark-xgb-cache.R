# benchmark_xgb_level(.cache = env): the domain list and the domain -> benchmark-group crosswalk are cached across calls (megb benchmarks every
# bootstrap replicate against the same census). The benchmarked values must be identical to the uncached call, for every call that reuses the cache.
test_that("benchmark_xgb_level with a cache is identical to the uncached call", {
  set.seed(7)
  D <- 30; per <- 12
  pop <- data.frame(dom = rep(sprintf("d%02d", 1:D), each = per), w = rpois(D * per, 40) + 1, stringsAsFactors = FALSE)
  pop$grp <- paste0("g", (as.integer(sub("d", "", pop$dom)) - 1) %/% 6 + 1)
  smp <- data.frame(dom = rep(sprintf("d%02d", 1:D), each = 3), y = runif(D * 3), sw = runif(D * 3, 1, 3), stringsAsFactors = FALSE)
  smp$grp <- pop$grp[match(smp$dom, pop$dom)]
  fw <- list(pop_data = pop, smp_data = smp, domains = "dom", smp_weights = "sw", pop_weights = "w")
  bm <- setNames(c(0.30, 0.45, 0.50, 0.25, 0.60), paste0("g", 1:5))
  for (type in c("ratio", "logit_raking")) {
    cache <- new.env(parent = emptyenv())
    for (rep in 1:4) {
      est <- runif(D, 0.1, 0.8); lab <- sample(sprintf("d%02d", 1:D))   # replicate-specific estimates in a shuffled order, as in the bootstrap
      args <- list(point_estim = list(ind = data.frame(Mean = est, Domain = lab, stringsAsFactors = FALSE)), framework = fw,
                   fixed = y ~ 1, benchmark = bm, benchmark_type = type, benchmark_level = "grp")
      a <- do.call(povmap:::benchmark_xgb_level, args)
      b <- do.call(povmap:::benchmark_xgb_level, c(args, list(.cache = cache)))
      expect_identical(a, b, info = paste(type, "call", rep))
    }
    expect_false(is.null(cache$crosswalk)); expect_false(is.null(cache$pop_dom))
  }
})
