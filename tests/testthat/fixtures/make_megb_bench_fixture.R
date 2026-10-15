# Writes megb_bench_identity_v061.rds: ind, var and CI of megb() with a numeric benchmark (ratio and logit_raking, no
# perturbation) on the helper-megb-bench.R problem, from the v0.6.1 build (povmap 1.0.3, 3ca1d01). Run with that build
# loaded:  Rscript make_megb_bench_fixture.R <path to library holding povmap 1.0.3>
args <- commandArgs(TRUE); if (length(args)) .libPaths(c(args[1], .libPaths()))
suppressMessages(library(povmap)); stopifnot(as.character(packageVersion("povmap")) == "1.0.3")
here <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "..", "helper-megb-bench.R"))
out <- list()
for (nm in names(megb_bench_cfgs)) for (mt in c("var", "mse")) {
  m <- megb_bench_run(megb_bench_cfgs[[nm]], mse_type = mt)
  out[[paste(nm, mt, sep = "|")]] <- list(ind = m$ind, var = m$var, CI = m$CI)
}
saveRDS(out, file.path(here, "megb_bench_identity_v061.rds"))
cat("wrote", length(out), "scenarios\n")
