source("psi_validation/R/stochastic_psi.R")
stopifnot(requireNamespace("data.table", quietly = TRUE))
stopifnot(requireNamespace("Rcpp", quietly = TRUE))
Rcpp::sourceCpp("psi_validation/R/ctmc_psi.cpp")
dir.create("psi_validation/data", FALSE, TRUE)
args <- commandArgs(trailingOnly = TRUE)
checkpoint_every <- 1000L
ce_arg <- grep("^--checkpoint-every=", args, value = TRUE)
if (length(ce_arg)) checkpoint_every <- as.integer(sub("^--checkpoint-every=", "", ce_arg))
psi_values <- c(0.25, 0.5, 0.75, 1)
base_attempts <- 3000L
expanded_attempts <- 10000L
batch_size <- expanded_attempts
small <- expand.grid(
  rho = c(.01, .02, .05, .10), theta = c(0, .5, 1),
  K = c(1000, 3000, 10000, 30000), R0 = 1 + c(.05, .075, .10, .15, .20, .30, .50, .75, 1, 1.5, 2, 3, 5)
)
small$scan <- "exact"
small <- small[order(small$rho, small$theta, small$K, small$R0), ]
cl <- parallel::makeCluster(8L)
on.exit(parallel::stopCluster(cl), add = TRUE)
parallel::clusterExport(cl, c("one_ctmc_psi", "one_adaptive_tau_psi"), envir = environment())
for (psi_idx in seq_along(psi_values)) {
  psi <- psi_values[psi_idx]
  label <- sub("\\.", "", as.character(psi))
  grid <- small
  grid$psi <- psi
  grid$point_id <- seq_len(nrow(grid))
  checkpoint <- sprintf("psi_validation/data/psi%s_scan_checkpoint.rds", label)
  outfile <- sprintf("psi_validation/data/psi%s_stochastic_results.csv", label)
  version <- sprintf("psi%s_dynamic_g_v1", label)
  new_state <- function() {
    list(version = version, grid = grid, counts = transform(grid,
      attempts = 0L,
      established = 0L, persistent = 0L, unresolved = 0L, target_attempts = base_attempts, batches = 0L, elapsed_seconds = 0
    ))
  }
  state <- if (file.exists(checkpoint)) readRDS(checkpoint) else new_state()
  if (!identical(state$grid, grid)) {
    old <- state
    state <- new_state()
    ncopy <- min(nrow(grid), nrow(old$counts))
    stopifnot(all(old$grid$scan[seq_len(ncopy)] == "exact"))
    state$counts[seq_len(ncopy), ] <- old$counts[seq_len(ncopy), names(state$counts)]
  }
  save_state <- function() {
    saveRDS(state, paste0(checkpoint, ".tmp"), compress = FALSE)
    file.copy(paste0(checkpoint, ".tmp"), checkpoint, TRUE)
    file.remove(paste0(checkpoint, ".tmp"))
    z <- data.table::as.data.table(state$counts)
    z[, `:=`(
      P_unconditional = persistent / attempts,
      uncond_low = NA_real_, uncond_high = NA_real_, P_conditional = NA_real_, cond_low = NA_real_, cond_high = NA_real_
    )]
    for (j in which(z$attempts > 0)) {
      u <- wilson(z$persistent[j], z$attempts[j])
      z[j, `:=`(uncond_low = u[1], uncond_high = u[2])]
      if (z$established[j] > 0) {
        cc <- wilson(z$persistent[j], z$established[j])
        z[j, `:=`(P_conditional = persistent / established, cond_low = cc[1], cond_high = cc[2])]
      }
    }
    data.table::fwrite(z, outfile)
  }
  batches_done <- 0L
  maybe_save <- function() {
    batches_done <<- batches_done + 1L
    if (batches_done %% checkpoint_every == 0) save_state()
  }
  for (i in seq_len(nrow(grid))) {
    repeat{
      x <- state$counts[i, ]
      if (x$attempts >= x$target_attempts) {
        if (x$attempts == base_attempts) {
          wu <- diff(wilson(x$persistent, x$attempts))
          wc <- if (x$established > 0) diff(wilson(x$persistent, x$established)) else Inf
          if (max(wu, wc) > .05) {
            state$counts$target_attempts[i] <- expanded_attempts
            maybe_save()
            next
          }
        }
        break
      }
      n <- min(batch_size, x$target_attempts - x$attempts)
      seed <- 1950000000L + psi_idx * 10000000L + i * 100L + x$batches + 1L
      t0 <- proc.time()[["elapsed"]]
      got <- if (x$scan == "exact") {
        simulate_counts_cpp(x$R0, x$rho, x$theta, psi, x$K, n, seed)
      } else {
        simulate_counts_batch_tau_psi(cl, x$R0, x$rho, x$theta, psi, x$K, n, seed)
      }
      state$counts$attempts[i] <- x$attempts + got["attempts"]
      state$counts$established[i] <- x$established + got["established"]
      state$counts$persistent[i] <- x$persistent + got["persistent"]
      state$counts$unresolved[i] <- x$unresolved + got["unresolved"]
      state$counts$batches[i] <- x$batches + 1L
      state$counts$elapsed_seconds[i] <- x$elapsed_seconds + proc.time()[["elapsed"]] - t0
      maybe_save()
      cat(sprintf("psi=%s %d/%d %s %d/%d\n", label, i, nrow(grid), x$scan, state$counts$attempts[i], state$counts$target_attempts[i]))
      if (got["unresolved"] > 0) stop("unresolved at psi=", label, " point ", i)
    }
  }
  save_state()
  cat(sprintf("psi=%s scan complete\n", label))
}
cat("all psi scans complete\n")
