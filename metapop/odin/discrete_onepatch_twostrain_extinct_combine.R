## Combine per-task outputs from discrete_onepatch_twostrain_extinct_run_array.R
library(dplyr)
library(tidyr)

files <- sort(list.files("outputs", pattern = "^twostrain_extinct_task_[0-9]+\\.rds$",
                         full.names = TRUE))
if (length(files) == 0) stop("no task output files found in outputs/")

res_df <- bind_rows(lapply(files, readRDS))
dd2    <- res_df |> pivot_longer(cols = c("I1", "I2"))
dd2_i1 <- dd2 |> dplyr::filter(name == "I1")

saveRDS(dd2_i1, "discrete_onepatch_twostrain_extinct.rds")
cat("combined", nrow(res_df), "rows from", length(files), "files\n")
