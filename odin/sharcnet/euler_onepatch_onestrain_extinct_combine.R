## Combine per-task outputs from euler_onepatch_onestrain_extinct_run_array.R
## Usage: Rscript euler_onepatch_onestrain_extinct_combine.R <combination>
## <combination> is one of: logistic_continuous, logistic_reedfrost,
##                           linear_continuous, linear_reedfrost,
##                           logistic_continuous_demoggrid
##                           (demoggrid: grid also expanded over r)
library(dplyr)

valid <- c("logistic_continuous", "logistic_reedfrost",
           "linear_continuous",   "linear_reedfrost",
           "logistic_continuous_demoggrid")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1 || !args[1] %in% valid)
  stop("Usage: Rscript euler_onepatch_onestrain_extinct_combine.R <combination>\n",
       "  combination must be one of: ", paste(valid, collapse = ", "))

combo <- args[1]
base  <- paste0("euler_onepatch_onestrain_extinct_", combo)

files <- sort(list.files("outputs",
                         pattern = paste0("^", base, "_task_[0-9]+\\.rds$"),
                         full.names = TRUE))
if (length(files) == 0)
  stop("no task output files found in outputs/ for combination: ", combo)

out <- bind_rows(lapply(files, readRDS))
outfile <- file.path("outputs", paste0(base, ".rds"))
saveRDS(out, outfile)
cat("combined", nrow(out), "rows from", length(files), "files ->", outfile, "\n")
