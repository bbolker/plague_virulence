## Combine per-task outputs from pairwise_invasion_factorial_run_array.R
## into pairwise_invasion_factorial[_mini].rds, a list with
##   run_summary    one row per run (task_id, parameters, replicate, outcomes)
##   combo_summary  one row per (K, alpha, seasonal_amp, R0_res, R0_inv)
##   thinned        occupancy series of the saved replicates
##   settings       fixed settings from pairwise_invasion_functions.R
## Missing task ids are reported (and listed in the output) but don't stop it.
library(dplyr)
library(here)
library(optparse)

opt <- parse_args(OptionParser(option_list = list(
  make_option(c("-m", "--mini"), action = "store_true", default = FALSE,
              help = "combine mini task outputs")
)))

source(here::here("fadeout", "seasonal", "pairwise_invasion_functions.R"))

base_fn <- "pairwise_factorial"
if (opt$mini) base_fn <- paste0(base_fn, "_mini")
n_tasks <- if (opt$mini) 32L else 1323L

cat("starting\n")
pattern <- sprintf("^%s_task_[0-9]+\\.rds$", base_fn)
files <- sort(list.files("outputs", pattern = pattern, full.names = TRUE))
if (length(files) == 0) stop("no task output files found in outputs/ matching ", pattern)
cat(length(files), "files found\n")

found <- as.integer(sub(".*_task_([0-9]+)\\.rds$", "\\1", files))
missing <- setdiff(seq_len(n_tasks), found)
if (length(missing)) {
  cat(length(missing), "missing task ids:", head(missing, 50),
      if (length(missing) > 50) "...", "\n")
}

res <- lapply(files, readRDS)
run_summary <- bind_rows(lapply(res, `[[`, "summary"))
thinned <- bind_rows(lapply(res, `[[`, "thinned"))
rm(res)

combo_summary <- summarise_pairs(run_summary, K, alpha, seasonal_amp,
                                 R0_res, R0_inv)

settings <- mget(c("n_patch", "resident_years", "invasion_years",
                   "summary_years", "I_invader", "THIN", "save_from_year"))

out_fn <- paste0(base_fn, ".rds")
cat("saving ...\n")
saveRDS(list(run_summary = run_summary, combo_summary = combo_summary,
             thinned = thinned, settings = settings, missing_tasks = missing),
        out_fn)
cat("combined", nrow(run_summary), "runs from", length(files), "files into",
    out_fn, "\n")
