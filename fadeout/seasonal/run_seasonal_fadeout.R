## Command-line runner for seasonal recurrent-fade-out diagnostics.
## Run from the repository root; see fadeout/seasonal/README.md.

library(here)
library(optparse)

source(here::here(
  "fadeout", "seasonal", "seasonal_fadeout_functions.R"
))

options <- list(
  make_option(
    "--grid", type = "character",
    default = here::here(
      "fadeout", "seasonal", "seasonal_fadeout_example_grid.csv"
    ),
    help = "Parameter-grid CSV"
  ),
  make_option(
    "--run-id", type = "character", default = "example",
    help = "Short output identifier [default %default]"
  ),
  make_option(
    "--output-root", type = "character",
    default = here::here("fadeout", "output", "seasonal_fadeout"),
    help = "Output root directory"
  ),
  make_option(
    "--dry-run", action = "store_true", default = FALSE,
    help = "Validate and print the grid without simulating"
  )
)
opt <- parse_args(OptionParser(option_list = options))

if (!file.exists(opt$grid)) stop("Grid file not found: ", opt$grid)
if (!grepl("^[A-Za-z0-9_.-]+$", opt$`run-id`)) {
  stop("run-id may contain only letters, numbers, dot, underscore, and hyphen")
}
grid <- validate_seasonal_grid(utils::read.csv(
  opt$grid, stringsAsFactors = FALSE, check.names = FALSE
))

cat("Validated parameter combinations: ", nrow(grid), "\n", sep = "")
cat("Total stochastic replicates: ", sum(grid$n_reps), "\n", sep = "")
print(grid)
if (opt$`dry-run`) quit(save = "no", status = 0)

output_dir <- file.path(opt$`output-root`, opt$`run-id`)
model_file <- here::here(
  "fadeout", "seasonal", "seasonal_model_metapop.R"
)
run_seasonal_fadeout_grid(grid, output_dir, model_file)
cat("Output: ", output_dir, "\n", sep = "")
