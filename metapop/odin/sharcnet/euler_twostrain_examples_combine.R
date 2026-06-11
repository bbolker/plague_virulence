## Combine per-task outputs from euler_twostrain_examples_run_array.R
## Run from the metapop/odin/sharcnet directory.
## Saves a list of 6 thinned tibbles (one per parameter set) matching the
## format of euler_twostrain_singlepatchintro_examples.R.
library(dplyr)

pars    <- read.csv("euler_twostrain_example_pars.csv")
pattern <- "^euler_twostrain_examples_task_[0-9]+\\.rds$"
files   <- sort(list.files("outputs", pattern = pattern, full.names = TRUE))
if (length(files) == 0) stop("no task output files found in outputs/ matching ", pattern)

cat(length(files), "files found\n")
results <- lapply(files, readRDS)
attr(results, "metadata") <- pars

saveRDS(results, "euler_twostrain_examples.rds")
cat("saved", length(results), "parameter sets\n")
