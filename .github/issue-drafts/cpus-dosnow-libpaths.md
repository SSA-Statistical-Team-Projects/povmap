# `xgb(cpus > 1)` fails with "object '.doSnowGlobals' not found" when the library is not on the workers' path

## Summary

`xgb()` with `cpus > 1` aborts during the bootstrap:

```
Error in checkForRemoteErrors(lapply(cl, recvResult)) :
  6 nodes produced errors; first error: object '.doSnowGlobals' not found
```

The cause is environmental rather than a logic error in povmap, but the failure
mode is opaque and povmap can either prevent it or diagnose it. Filing so the
behaviour is recorded.

## Root cause

`R/xgb.R:570-571`:

```r
cl <- parallel::makeCluster(cpus)
doSNOW::registerDoSNOW(cl)
```

`parallel::makeCluster()` starts **PSOCK** workers, which are fresh R sessions.
They inherit environment variables but **not** a `.libPaths()` set at runtime in
the parent. If povmap and its dependencies live in a library added at runtime
(`.libPaths(c(lib, .libPaths()))`, common under `renv`, a shared team library,
or a non-default `R_LIBS` layout), the workers cannot see `doSNOW`.
`registerDoSNOW()` then fails to initialise `.doSnowGlobals` on them, and the
first `%dopar%` iteration reports the missing object rather than the missing
package.

Demonstrated directly:

```r
.libPaths(c("E:/David/R/lib-4.5.3", .libPaths()))
cl <- parallel::makeCluster(2)
parallel::clusterEvalQ(cl, .libPaths())
#> [1] "G:/David/R/R-4.5.3/library"          # runtime path absent
parallel::clusterEvalQ(cl, "doSNOW" %in% rownames(installed.packages()))
#> [1] FALSE FALSE

Sys.setenv(R_LIBS_USER = "E:/David/R/lib-4.5.3")   # env vars ARE inherited
cl2 <- parallel::makeCluster(2)
parallel::clusterEvalQ(cl2, "doSNOW" %in% rownames(installed.packages()))
#> [1] TRUE TRUE
doSNOW::registerDoSNOW(cl2)                        # succeeds
```

With `R_LIBS_USER` set, `xgb(cpus = 6, B = 100)` runs to completion.

## Not platform-specific in principle, but Windows makes it certain

PSOCK is the only cluster type on Windows, so the problem always applies there.
On Linux/macOS the default is still PSOCK for `parallel::makeCluster()`, so the
same failure occurs; a `type = "FORK"` cluster would inherit the parent's
`.libPaths()` and mask it.

## `megb()` is unaffected

`megb()` and `megb_mse()` contain no cluster construction — no `makeCluster`,
`registerDoSNOW`, `registerDoParallel` or `%dopar%`. The bug is confined to
`xgb()`'s bootstrap.

## Suggested fixes, in order of preference

1. Propagate the parent's library path to the workers immediately after the
   cluster is created, before `registerDoSNOW()`:
   ```r
   cl <- parallel::makeCluster(cpus)
   parallel::clusterCall(cl, function(p) .libPaths(p), .libPaths())
   doSNOW::registerDoSNOW(cl)
   ```
   Two lines, no new dependency, fixes every runtime-library layout.

2. Failing that, check reachability and fail with a useful message:
   ```r
   if (!all(unlist(parallel::clusterEvalQ(cl, requireNamespace("doSNOW", quietly = TRUE))))) {
     parallel::stopCluster(cl)
     stop("cpus > 1 requires doSNOW on the workers' .libPaths(); ",
          "set R_LIBS_USER or use cpus = 1.")
   }
   ```

## Separate observation: `cpus > 1` can be much slower

Measured on this workload (428,726 population sub-areas, 84,962 sample rows,
1,122 domains, `nrounds = 100`):

| B | cpus | wall |
|---|------|------|
| 100 | 1 | 17.5 s |
| 100 | 6 | 56.0 s |

The bootstrap itself is ~8 s of the 17.5 s; the rest is data preparation, which
is not parallelised. Six worker sessions each loading povmap and its ~37
imports costs more than the work they receive. Worth a note in `?xgb` that
`cpus > 1` pays off only when `B` and the per-iteration cost are large enough to
amortise worker startup.
