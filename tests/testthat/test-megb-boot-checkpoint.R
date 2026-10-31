# Opt-in progress log and checkpoint / resume of the leaves_only bootstrap (R/megb_boot_checkpoint.R). The replicate draws no random numbers, so an
# interrupted run that is resumed from its checkpoint must give the same numbers as an uninterrupted run (B = 50 here). Needs the installed package
# (the parallel runs load it on the workers). The "stop_after" option simulates an interruption after that many replicates; the last test kills a real process.
skip_ck <- function() {
  skip_on_cran(); skip_if_not_installed("xgboost"); skip_if_not_installed("collapse"); skip_if_not_installed("jsonlite"); skip_if_not_installed("rlang")
  skip_if_not("povmap" %in% rownames(utils::installed.packages()))
  old <- Sys.getenv("OMP_NUM_THREADS"); Sys.setenv(OMP_NUM_THREADS = "1")
  function() Sys.setenv(OMP_NUM_THREADS = old)
}
chunk_files <- function(d) list.files(d, pattern = "^chunk_.*[.]rds$")

for (cores in c(0, 2)) {
  test_that(paste0("an interrupted and resumed bootstrap is identical to an uninterrupted one (", if (cores) "parallel" else "sequential", ")"), {
    r <- skip_ck(); on.exit(r(), add = TRUE)
    d <- ck_data(); dir <- tempfile("ck"); on.exit(unlink(dir, recursive = TRUE), add = TRUE)
    ref <- suppressMessages(ck_fit(cores, d = d))
    o <- options(povmap.megb_boot.stop_after = 20)
    expect_error(suppressMessages(ck_fit(cores, d = d, boot_checkpoint_dir = dir, boot_chunk = 10)), "simulated interruption after 20")
    options(o)
    expect_equal(length(chunk_files(dir)), 2L)                                    # 20 of 50 replicates saved
    msg <- capture_messages(res <- ck_fit(cores, d = d, boot_checkpoint_dir = dir, boot_chunk = 10))
    expect_true(any(grepl("20 of 50 replicates already done, 30 to run", msg)))
    expect_equal(length(chunk_files(dir)), 5L)
    ck_same(ref, res, paste("resumed", cores))
    again <- suppressMessages(ck_fit(cores, d = d, boot_checkpoint_dir = dir, boot_chunk = 10))   # a finished checkpoint: nothing left to run
    ck_same(ref, again, paste("finished checkpoint", cores))
  })
}

test_that("a run interrupted sequentially can be resumed on a cluster, and with another chunk size", {
  r <- skip_ck(); on.exit(r(), add = TRUE)
  d <- ck_data(); dir <- tempfile("ck"); on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  ref <- suppressMessages(ck_fit(0, d = d))
  o <- options(povmap.megb_boot.stop_after = 20)
  expect_error(suppressMessages(ck_fit(0, d = d, boot_checkpoint_dir = dir, boot_chunk = 10)), "simulated interruption")
  options(o)
  res <- suppressMessages(ck_fit(2, d = d, boot_checkpoint_dir = dir, boot_chunk = 7))
  ck_same(ref, res, "sequential then parallel")
})

test_that("the defaults are off: no files, no progress lines, and the same numbers", {
  r <- skip_ck(); on.exit(r(), add = TRUE)
  d <- ck_data()
  expect_null(formals(povmap::megb)$boot_checkpoint_dir); expect_false(formals(povmap::megb)$boot_progress)
  expect_null(povmap:::.megb_boot_ckpt(NULL, 50, NULL, NULL, FALSE, NULL, 0))
  msg <- capture_messages(a <- ck_fit(0, d = d)); expect_false(any(grepl("estimated remaining", msg)))
  pf <- tempfile(fileext = ".log"); on.exit(unlink(pf), add = TRUE)
  msg2 <- capture_messages(b <- ck_fit(0, d = d, boot_progress = TRUE, boot_progress_file = pf, boot_chunk = 10))   # progress alone changes nothing else
  ck_same(a, b, "progress only")
  expect_true(any(grepl("50 of 50 replicates \\(100.0%\\).*estimated remaining 0:00:00", msg2)))
  lines <- readLines(pf); expect_true(length(lines) >= 6L); expect_true(all(grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2} ", lines)))
  expect_true(any(grepl("20 of 50 replicates \\(40.0%\\), 20 run in this session; elapsed [0-9:]+, estimated remaining [0-9:]+", lines)))
})

test_that("progress is logged from the main process on a cluster", {
  r <- skip_ck(); on.exit(r(), add = TRUE)
  msg <- capture_messages(ck_fit(2, d = ck_data(), boot_progress = TRUE, boot_chunk = 25))
  expect_true(any(grepl("25 of 50 replicates", msg))); expect_true(any(grepl("50 of 50 replicates", msg)))
})

test_that("a checkpoint written by a different call or different data is refused", {
  r <- skip_ck(); on.exit(r(), add = TRUE)
  d <- ck_data(); dir <- tempfile("ck"); on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  suppressMessages(ck_fit(0, d = d, boot_checkpoint_dir = dir, boot_chunk = 25))
  n <- length(chunk_files(dir)); expect_equal(n, 2L)
  expect_error(suppressMessages(ck_fit(0, d = d, seed = 12, boot_checkpoint_dir = dir)), "refusing to resume")        # other seed: other draws
  expect_error(suppressMessages(ck_fit(0, d = d, B = 60, boot_checkpoint_dir = dir)), "refusing to resume.*B: saved 50, now 60")
  d2 <- d; d2$smp$x1[3] <- d2$smp$x1[3] + 1e-6                                   # one covariate value in the sample
  expect_error(suppressMessages(ck_fit(0, d = d2, boot_checkpoint_dir = dir)), "refusing to resume.*(smp_data|X_smp_mat)")
  d3 <- d; d3$pop$x2[5] <- d3$pop$x2[5] + 1e-6                                   # one census covariate value
  expect_error(suppressMessages(ck_fit(0, d = d3, boot_checkpoint_dir = dir)), "refusing to resume.*X_pop_mat")
  d4 <- d; d4$bm <- d4$bm * 1.001                                                # the benchmark targets
  expect_error(suppressMessages(ck_fit(0, d = d4, boot_checkpoint_dir = dir)), "refusing to resume")
  expect_error(suppressMessages(ck_fit(0, d = d, weightedBS = FALSE, boot_checkpoint_dir = dir)), "refusing to resume")
  expect_error(suppressMessages(ck_fit(0, d = d, em_maxit = 30L, boot_checkpoint_dir = dir)), "refusing to resume")    # a setting that changes the fit
  expect_equal(length(chunk_files(dir)), n)                                      # a refused call writes nothing
  .real_pkg_id <- povmap:::.megb_pkg_id
  # another povmap / xgboost / lme4 build (mocked in the main process; the fit is sequential)
  testthat::local_mocked_bindings(.megb_pkg_id = function() { x <- .real_pkg_id(); x$povmap_version <- "9.9.9"; x }, .package = "povmap")
  expect_error(suppressMessages(ck_fit(0, d = d, boot_checkpoint_dir = dir)), "refusing to resume.*pkg")
  testthat::local_mocked_bindings(.megb_pkg_id = function() { x <- .real_pkg_id(); x$xgboost_version <- "0.0.1"; x }, .package = "povmap")
  expect_error(suppressMessages(ck_fit(0, d = d, boot_checkpoint_dir = dir)), "refusing to resume.*pkg")
  testthat::local_mocked_bindings(.megb_pkg_id = .real_pkg_id, .package = "povmap")
  # the same call with another chunk size is accepted
  expect_no_error(suppressMessages(ck_fit(0, d = d, boot_checkpoint_dir = dir, boot_chunk = 5)))
  # damaged directories
  f <- file.path(dir, chunk_files(dir)[1]); x <- readRDS(f); x$hash <- "other"; saveRDS(x, f)
  expect_error(suppressMessages(ck_fit(0, d = d, boot_checkpoint_dir = dir)), "does not belong to this call")
  file.remove(file.path(dir, "manifest.rds"))
  expect_error(suppressMessages(ck_fit(0, d = d, boot_checkpoint_dir = dir)), "no manifest.rds")
})

test_that("the options are validated", {
  r <- skip_ck(); on.exit(r(), add = TRUE)
  d <- ck_data()
  expect_error(ck_fit(0, d = d, boot_checkpoint_dir = c("a", "b")), "single directory path")
  expect_error(suppressMessages(ck_fit(0, d = d, boot_checkpoint_dir = tempfile(), boot_chunk = 0)), "boot_chunk")
  expect_error(suppressMessages(povmap::megb(fixed = y ~ x1 + x2, smp_data = d$smp, pop_data = d$pop, domains = "dom", mse = FALSE, boot_progress = TRUE)), "need mse = TRUE")
})

test_that("a real kill after about 20 replicates, then a resume, gives the uninterrupted result (sequential and parallel)", {
  r <- skip_ck(); on.exit(r(), add = TRUE)
  skip_if_not_installed("processx")
  helper <- normalizePath(testthat::test_path("helper-megb-boot-checkpoint.R"), winslash = "/")
  rs <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
  libs <- paste0('"', normalizePath(.libPaths(), winslash = "/"), '"', collapse = ", ")
  d <- ck_data(); ref <- suppressMessages(ck_fit(0, d = d))
  for (cores in c(0, 2)) {
    dir <- tempfile("ck"); script <- tempfile(fileext = ".R")
    writeLines(c(sprintf('.libPaths(c(%s))', libs), sprintf('source("%s")', helper), 'Sys.setenv(OMP_NUM_THREADS = "1")',
                 sprintf('x <- ck_fit(%d, boot_checkpoint_dir = "%s", boot_chunk = 5)', as.integer(cores), gsub("\\\\", "/", dir))), script)
    p <- processx::process$new(rs, c("--vanilla", script), stdout = "|", stderr = "|")
    t0 <- Sys.time()
    while (p$is_alive() && length(chunk_files(dir)) < 4L && difftime(Sys.time(), t0, units = "secs") < 120) Sys.sleep(0.05)   # 4 chunks of 5 = 20 replicates
    killed <- p$is_alive(); p$kill_tree(); p$wait()
    expect_true(killed, info = "the child finished before it was killed; slow the fit down")
    n <- length(chunk_files(dir)); expect_true(n >= 4L && n < 10L, info = paste(n, "chunks after the kill"))
    res <- suppressMessages(ck_fit(cores, d = d, boot_checkpoint_dir = dir, boot_chunk = 5))
    ck_same(ref, res, paste("real kill, cores", cores))
    unlink(c(dir, script), recursive = TRUE)
  }
})
