## Combine per-task outputs from euler_onestrain_run_array.R
library(dplyr)
library(optparse)

opt <- parse_args(OptionParser(option_list = list(
  make_option(c("-m", "--mini"),   action = "store_true", default = FALSE,
              help = "combine mini task outputs"),
  make_option(c("-b", "--batch2"), action = "store_true", default = FALSE,
              help = "combine batch2 task outputs")
)))

base_fn <- "euler_onestrain_intro"
if (opt$mini)   base_fn <- paste0(base_fn, "_mini")
if (opt$batch2) base_fn <- paste0(base_fn, "_batch2")

pattern <- sprintf("^%s_task_[0-9]+\\.rds$", base_fn)
files <- sort(list.files("outputs", pattern = pattern, full.names = TRUE))
if (length(files) == 0) stop("no task output files found in outputs/ matching ", pattern)

out <- bind_rows(lapply(files, readRDS))
saveRDS(out, paste0(base_fn, ".rds"))
cat("combined", nrow(out), "rows from", length(files), "files\n")
