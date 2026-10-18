# PSOCK workers started by povmap (cpus > 1) load povmap from their own library path. Pin it to the library the package under test
# was loaded from, so a worker can never pick up another installed povmap (an older one lacks internals such as .boot_domain_wmean).
# R_LIBS is inherited by the worker processes. Skipped when the package is not loaded from an installed library (load_all).
local({
  pkg <- find.package("povmap")
  if (file.exists(file.path(pkg, "Meta", "package.rds"))) {
    old <- Sys.getenv("R_LIBS", unset = NA)
    libs <- unique(c(normalizePath(dirname(pkg), winslash = "/", mustWork = FALSE), .libPaths()))
    Sys.setenv(R_LIBS = paste(libs, collapse = .Platform$path.sep))
    if (requireNamespace("withr", quietly = TRUE)) {
      withr::defer(if (is.na(old)) Sys.unsetenv("R_LIBS") else Sys.setenv(R_LIBS = old), envir = testthat::teardown_env())
    }
  }
})
