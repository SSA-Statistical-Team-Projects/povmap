# Progress log and checkpoint / resume for the leaves_only bootstrap of mse_megb() (opt-in; all defaults off, existing calls unchanged).
#
# The leaves_only replicate draws no random numbers (y_star_smp, the bootstrap truth and the benchmark targets are made before the loop), so
# replicate i is a function of the fit and of column i of y_star_smp alone. Finished replicates can therefore be saved and, on a rerun with the
# same inputs, loaded instead of recomputed, and the collated result is identical to an uninterrupted run (tests/testthat/test-megb-boot-checkpoint.R).
#
# boot_checkpoint_dir: directory. manifest.rds holds a fingerprint of the call and the data; chunk_<first>_<n>.rds holds the replicates finished in
#   one chunk (written to a .tmp file and renamed, so a kill during a write leaves no partial chunk). A rerun with a different fingerprint stops.
# boot_chunk: replicates per chunk (checkpoint granularity and progress interval). Default 25 sequentially and max(25, 4 * workers) on a cluster:
#   a cluster chunk is one parLapplyLB() call, so workers wait for the slowest replicate at each chunk end; smaller chunks cost idle time.
# boot_progress / boot_progress_file: a line per chunk from the main process (a worker's messages are lost under PSOCK) with the replicates
#   done, elapsed time and the estimated time remaining, to the console and / or appended to the file.

.megb_boot_fmt_time <- function(s) {
  if (!is.finite(s)) return("NA")
  s <- round(s); sprintf("%d:%02d:%02d", s %/% 3600, (s %% 3600) %/% 60, s %% 60)
}

# What goes into the fingerprint is made hashable: boosters are external pointers, so they enter as their saved bytes; a closure enters as its
# source plus the plain data objects in its enclosing frame (a benchmark function carries its targets there).
.megb_hashable <- function(x) {
  if (inherits(x, "xgb.Booster")) return(list(booster = xgboost::xgb.save.raw(x, raw_format = "ubj")))
  if (is.function(x)) {
    e <- environment(x); keep <- list()
    if (!is.null(e) && !identical(e, globalenv()) && !isNamespace(e) && !identical(e, baseenv())) {
      nm <- sort(ls(e, all.names = TRUE))
      keep <- mget(nm, envir = e)
      keep <- keep[vapply(keep, function(v) is.atomic(v) || is.data.frame(v), NA)]
    }
    return(list(src = deparse(x), data = keep, attrs = attributes(x)[setdiff(names(attributes(x)), "srcref")]))
  }
  if (is.list(x) && !is.data.frame(x)) { at <- attributes(x); x <- lapply(x, .megb_hashable); attributes(x) <- at[setdiff(names(at), "srcref")]; return(x) }
  x
}

.megb_boot_fingerprint <- function(parts) {
  h <- vapply(parts, function(p) rlang::hash(.megb_hashable(p)), "")
  list(hash = rlang::hash(h), parts = h)
}

.megb_boot_rds_write <- function(obj, path) {
  tmp <- paste0(path, ".tmp"); saveRDS(obj, tmp); if (!file.rename(tmp, path)) stop("could not write ", path, call. = FALSE)
}

# Returns NULL when neither progress nor checkpointing is asked for. Otherwise a list of functions sharing one state:
#   $done (indices already finished and loaded), $todo (indices still to run), $loaded (list of B, filled at $done),
#   $save(idx, res)  records a finished chunk (checkpoint file when a directory is set, then the progress line),
#   $chunk (size), $begin(i, boots) / $finish(boots) for the sequential loop.
.megb_boot_ckpt <- function(dir, B, parts, chunk, progress, progress_file, cores) {
  if (is.null(dir) && !isTRUE(progress) && is.null(progress_file)) return(NULL)
  if (!is.null(chunk)) {
    if (length(chunk) != 1L || !is.numeric(chunk) || !is.finite(chunk) || chunk < 1) stop("boot_chunk must be a single number >= 1.", call. = FALSE)
    chunk <- as.integer(chunk)
  } else chunk <- if (!is.null(cores) && cores > 1) max(25L, 4L * as.integer(cores)) else 25L
  st <- new.env(parent = emptyenv())
  st$loaded <- vector("list", B); st$done <- integer(0); st$t0 <- Sys.time(); st$ran <- 0L
  st$pend <- integer(0); st$cur <- NULL; st$stop_after <- getOption("povmap.megb_boot.stop_after")   # st$stop_after: test hook (simulated interruption)

  if (!is.null(dir)) {
    fp <- .megb_boot_fingerprint(parts)
    dir.create(dir, showWarnings = FALSE, recursive = TRUE)
    if (!dir.exists(dir)) stop("boot_checkpoint_dir ", dir, " could not be created.", call. = FALSE)
    mf <- file.path(dir, "manifest.rds")
    man <- list(format = 1L, B = B, hash = fp$hash, parts = fp$parts, created = format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
    stale <- list.files(dir, pattern = "[.]tmp$", full.names = TRUE); if (length(stale)) unlink(stale)   # a write that a kill cut short
    files <- list.files(dir, pattern = "^chunk_[0-9]+_[0-9]+[.]rds$", full.names = TRUE)
    if (file.exists(mf)) {
      old <- tryCatch(readRDS(mf), error = function(e) NULL)
      if (is.null(old) || !identical(old$format, 1L)) stop("boot_checkpoint_dir ", dir, ": manifest.rds is unreadable or of another format; refusing to resume. Use an empty directory.", call. = FALSE)
      if (!identical(old$B, B) || !identical(old$hash, fp$hash)) {
        o <- unname(old$parts[names(fp$parts)]); diff <- names(fp$parts)[is.na(o) | o != fp$parts]
        diff <- unique(c(diff, setdiff(names(old$parts), names(fp$parts))))
        stop("boot_checkpoint_dir ", dir, " was written by a different call (created ", old$created, "); refusing to resume.\n  B: saved ", old$B, ", now ", B,
             "\n  parts that differ: ", if (length(diff)) paste(diff, collapse = ", ") else "(none by name)",
             "\n  Use an empty directory, or delete this one to start again.", call. = FALSE)
      }
    } else {
      if (length(files)) stop("boot_checkpoint_dir ", dir, " holds replicate files but no manifest.rds; refusing to use it. Use an empty directory.", call. = FALSE)
      .megb_boot_rds_write(man, mf)
    }
    for (f in files) {
      ck <- tryCatch(readRDS(f), error = function(e) NULL)
      if (is.null(ck) || !identical(ck$hash, fp$hash) || !is.integer(ck$idx) || length(ck$idx) != length(ck$res) || any(ck$idx < 1L | ck$idx > B) || any(ck$idx %in% st$done))
        stop("boot_checkpoint_dir ", dir, ": ", basename(f), " is unreadable, overlaps another chunk or does not belong to this call; refusing to resume.", call. = FALSE)
      st$loaded[ck$idx] <- ck$res; st$done <- c(st$done, ck$idx)
    }
    st$done <- sort(st$done)
    if (length(st$done)) message("megb bootstrap: resuming from ", dir, ": ", length(st$done), " of ", B, " replicates already done, ", B - length(st$done), " to run")
  }
  st$todo <- setdiff(seq_len(B), st$done); st$chunk <- chunk; st$hash <- if (!is.null(dir)) fp$hash

  say <- function(txt) {
    if (isTRUE(progress)) message(txt)
    if (!is.null(progress_file)) try(cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), txt, "\n", file = progress_file, append = TRUE), silent = TRUE)
  }
  if (!is.null(progress_file) || isTRUE(progress)) say(sprintf("megb bootstrap: %d replicates, %d to run, chunks of %d%s", B, length(st$todo), chunk, if (length(st$done)) sprintf(" (%d resumed from checkpoint)", length(st$done)) else ""))
  st$save <- function(idx, res) {
    idx <- as.integer(idx)
    if (!is.null(dir)) .megb_boot_rds_write(list(hash = st$hash, idx = idx, res = res), file.path(dir, sprintf("chunk_%06d_%d.rds", idx[1], length(idx))))
    st$loaded[idx] <- res; st$done <- c(st$done, idx); st$ran <- st$ran + length(idx)
    el <- as.numeric(difftime(Sys.time(), st$t0, units = "secs")); left <- B - length(st$done)
    say(sprintf("megb bootstrap: %d of %d replicates (%.1f%%), %d run in this session; elapsed %s, estimated remaining %s", length(st$done), B, 100 * length(st$done) / B, st$ran,
                .megb_boot_fmt_time(el), .megb_boot_fmt_time(if (st$ran > 0) el / st$ran * left else NA_real_)))
    if (!is.null(st$stop_after) && st$ran >= st$stop_after && left > 0) stop("simulated interruption after ", st$ran, " replicates", call. = FALSE)
    invisible(NULL)
  }
  # sequential loop: begin(i) at the top of iteration i, finish() after the loop; the replicates before i are complete when begin(i) is called
  st$begin <- function(i, boots) {
    if (!is.null(st$cur)) st$pend <- c(st$pend, st$cur)
    st$cur <- i
    if (length(st$pend) >= st$chunk) { st$save(st$pend, boots[st$pend]); st$pend <- integer(0) }
  }
  st$finish <- function(boots) {
    st$pend <- c(st$pend, st$cur); st$cur <- NULL
    if (length(st$pend)) { st$save(st$pend, boots[st$pend]); st$pend <- integer(0) }
  }
  st
}
