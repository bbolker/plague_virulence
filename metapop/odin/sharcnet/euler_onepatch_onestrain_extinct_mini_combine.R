## Combine per-task outputs from euler_onepatch_onestrain_extinct_run_array.R
library(dplyr)

files <- sort(list.files("outputs", pattern = "^euler_onepatch_onestrain_extinct_mini_task_[0-9]+\\.rds$",
                         full.names = TRUE))
if (length(files) == 0) stop("no task output files found in outputs/")

out <- bind_rows(lapply(files, readRDS))
saveRDS(out, "outputs/euler_onepatch_onestrain_extinct_mini.rds")
cat("combined", nrow(out), "rows from", length(files), "files\n")
