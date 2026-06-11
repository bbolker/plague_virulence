library(plagueMetapop)
library(dplyr)
library(arrow)
## download big file:
##  scp nibi.sharcnet.ca:~/project/bolker/plague_virulence/metapop/odin/outputs/euler_twostrain_examples_task_000003.rds .

fn <- "metapop/odin/outputs/euler_twostrain_examples_task_000003.rds"
x <- readRDS(here::here(fn))

dt <- 0.1
strain2_delay <- 50

## delay
xmin <- x |> select(step:run) |>
    filter(step > strain2_delay) |>
    mutate(across(run, as.integer))
    
saveRDS(xmin, "tmp.rds")

write_parquet(xmin, "tmp.pqt") ## much smaller, FWIW
xmin <- read_parquet("tmp.pqt", as.data.frame = FALSE)

## not too much help (due to compression?)
file.size("tmp.rds")/file.size(fn)
file.size("tmp.pqt")/file.size(fn)

x <- xmin |> filter(run==1L)
x2 <- x |> filter(step == min(step))

sumval <- xmin |> summarise(patches = sum(value>0),
                mean = mean(value[value>0]),
                .by = c(state, step, run))
