# Real-data reproduction of xgb() estimation (Nigeria 2023, poverty, ward level). Skipped unless
# POVMAP_NGA_REPRO_ARGS (the xgb() argument list of the pov ward tip-control run: baseline features,
# the saved baseline tuning passed directly, no retune) and POVMAP_NGA_REPRO_EXPECTED (that run's
# saved result, povmap 315e705 + xgboost 3.1.2.1) point to files. Takes about 40 minutes on 2 cores.
# Guards that changes elsewhere in the package (e.g. xgb_tune) leave xgb() estimation unchanged.

test_that("xgb() reproduces the Nigeria pov ward tip-control estimates exactly", {
  args_file <- Sys.getenv("POVMAP_NGA_REPRO_ARGS")
  exp_file  <- Sys.getenv("POVMAP_NGA_REPRO_EXPECTED")
  skip_if(!nzchar(args_file) || !nzchar(exp_file) || !file.exists(args_file) || !file.exists(exp_file),
          "Nigeria reproduction inputs not available")
  args <- readRDS(args_file)
  expected <- readRDS(exp_file)
  got <- suppressMessages(do.call(xgb, args))
  expect_identical(got$ind, expected$ind)
  expect_identical(got$var, expected$var)
  expect_identical(got$CI, expected$CI)
})
