library(plagueMetapop)
library(future)
library(ggplot2); theme_set(theme_bw())

nsim <- 100L

plan(multicore(workers = 10L))
set.seed(101)

args <- list(beta_vec = c(3, 0),
             K        = 1e6,
             r        = 0.125,
             n_patch  = 1,
             nt       = 500,
             alpha    = 0,
             I_init    = c(10, 0),
             stop_cond = NULL,
             platform = "odin")
rfun <- function(...) {
  set.seed(101)
  args <- modifyList(args, list(...))
  do.call(discrete_run, args)
}

list.files(path = system.file("odin", package = "plagueMetapop"))
run1 <- rfun(nsim = 100)
run2 <- rfun(def_file = "euler_odin_def.R", nsim = 100)
run3 <- rfun(def_file = "euler_odin_def.R", nsim = 100,
             reedfrost = 1)
run4 <- rfun(def_file = "ode_odin_def.R")
run5 <- rfun(def_file = "euler_odin_def.R", nsim = 100,
             K = 1e9)


## test ...
stopifnot(all.equal(run1, run3))
## don't need run3 (det/reedfrost are redundant)

run_list <- list(reedfrost = run1, euler = run2, euler_large = run5,
                 ode = run4)

runx <- purrr::map_dfr(run_list,
                       ~ dplyr::bind_rows(., .id = "run"),
                       .id = "model") |>
  dplyr::mutate(across(model, ~factor(., levels = names(run_list))),
                value = ifelse(grepl("large", model),
                               value/1000, value))
  
dplyr::filter(runx, state == "I1") |>
  ggplot(aes(step, value, colour = model, linewidth = model,
             alpha = model)) +
  geom_line(aes(group = interaction(run, model))) +
  scale_linewidth_manual(values = c(1, 1, 1, 2)) +
  scale_alpha_manual(values = c(0.4, 0.4, 0.4, 1)) +
  scale_y_log10() +
  scale_colour_brewer(palette="Dark2")

## ... discrepancy between ode and euler_large ... ???

